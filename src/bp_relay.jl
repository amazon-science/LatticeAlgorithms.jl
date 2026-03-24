# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
  
# Licensed under the Apache License, Version 2.0 (the "License").
# You may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

using Random
using SparseArrays

function _relay_prior_llr_vector(p::Real, num_bits::Int)
    if !(0 < p < 0.5)
        error("The bit-flip probability p has to satisfy 0 < p < 0.5.")
    end
    return fill(log((1 - p) / p), num_bits)
end

function _relay_prior_llr_vector(p::AbstractVector{<:Real}, num_bits::Int)
    if length(p) != num_bits
        error("The probability vector length has to equal the number of columns of H.")
    end
    if any(x -> !(0 < x < 0.5), p)
        error("All bit-flip probabilities have to satisfy 0 < p_j < 0.5.")
    end
    return log.((1 .- p) ./ p)
end

_relay_weight(error::AbstractVector{<:Integer}, llr_prior::AbstractVector{<:Real}) =
    sum(Int64.(error) .* llr_prior)

function _relay_sample_gamma_vector(
    rng::AbstractRNG,
    num_bits::Int,
    leg_idx::Int;
    gamma_schedule=nothing,
    first_leg_gamma::Real=0.35,
    gamma_center::Real=0.3655,
    gamma_width::Real=1.239,
)
    if !isnothing(gamma_schedule)
        if gamma_schedule isa AbstractMatrix
            if size(gamma_schedule, 2) != num_bits
                error("gamma_schedule has to have num_bits columns.")
            end
            if !(1 <= leg_idx <= size(gamma_schedule, 1))
                error("gamma_schedule does not contain the requested leg index.")
            end
            return Float64.(vec(gamma_schedule[leg_idx, :]))
        elseif gamma_schedule isa AbstractVector
            if !(1 <= leg_idx <= length(gamma_schedule))
                error("gamma_schedule does not contain the requested leg index.")
            end
            gamma_vec = gamma_schedule[leg_idx]
            if length(gamma_vec) != num_bits
                error("Each gamma vector in gamma_schedule has to have length num_bits.")
            end
            return Float64.(collect(gamma_vec))
        else
            error("gamma_schedule has to be a matrix or a vector of vectors.")
        end
    end

    if leg_idx == 1
        return fill(Float64(first_leg_gamma), num_bits)
    end

    half_width = gamma_width / 2
    lower = gamma_center - half_width
    upper = gamma_center + half_width
    return lower .+ (upper - lower) .* rand(rng, num_bits)
end

function _relay_build_graph(H::SparseMatrixCSC{Int64, Int64})
    num_checks, num_bits = size(H)

    check_to_bits = [Int64[] for _ in 1:num_checks]
    bit_to_checks = [Int64[] for _ in 1:num_bits]
    check_bit_pos = [Dict{Int64, Int64}() for _ in 1:num_checks]
    bit_check_pos = [Dict{Int64, Int64}() for _ in 1:num_bits]

    rows, cols, _ = findnz(H)
    for (check_idx, bit_idx) in zip(rows, cols)
        push!(check_to_bits[check_idx], bit_idx)
        check_bit_pos[check_idx][bit_idx] = length(check_to_bits[check_idx])

        push!(bit_to_checks[bit_idx], check_idx)
        bit_check_pos[bit_idx][check_idx] = length(bit_to_checks[bit_idx])
    end

    return check_to_bits, bit_to_checks, check_bit_pos, bit_check_pos
end

"""
    dmem_bp_decode(H::AbstractMatrix{<:Integer}, s::AbstractVector{<:Integer}, p;
                   max_iter::Int=size(H, 2),
                   initial_marginals=nothing,
                   gamma=0.0)

Decode the binary syndrome equation `H * e = s (mod 2)` using the
Disordered Memory Belief Propagation (DMem-BP) update rules introduced by
Müller et al. (IBM Quantum) in their Relay-BP work.

This implementation follows the min-sum check update together with the
memory-biased variable update
`Λ_j(t) = (1 - γ_j) Λ_j(0) + γ_j M_j(t - 1)`.

Arguments:
- `H`: binary parity-check matrix
- `s`: binary syndrome vector
- `p`: either a scalar bit-flip rate or a vector of per-bit rates

Keyword arguments:
- `max_iter`: maximum number of BP iterations
- `initial_marginals`: the initial marginals `M(0)`; by default these are the
  channel log-likelihood ratios
- `gamma`: either a scalar memory strength or a vector of per-bit memory
  strengths

The return value is a named tuple with fields
- `converged`: whether a syndrome-matching solution was found
- `error`: the final hard decision
- `marginals`: the final marginals
- `bias`: the final bias vector `Λ(t)`
- `iterations`: the number of iterations performed
- `syndrome_weight`: the Hamming weight of `H * error + s (mod 2)`
- `weight`: the weighted error cost `sum(error_j * log((1-p_j)/p_j))`
- `gamma`: the memory strengths used in this leg
"""
function dmem_bp_decode(
    H::AbstractMatrix{<:Integer},
    s::AbstractVector{<:Integer},
    p;
    max_iter::Int=size(H, 2),
    initial_marginals=nothing,
    gamma=0.0,
)
    H = sparse(mod.(Int64.(H), 2))
    s = mod.(Int64.(s), 2)

    num_checks, num_bits = size(H)
    if length(s) != num_checks
        error("The syndrome length has to equal the number of rows of H.")
    end
    if max_iter < 1
        error("max_iter has to be a positive integer.")
    end

    llr_prior = _relay_prior_llr_vector(p, num_bits)

    if isnothing(initial_marginals)
        marginals_prev = copy(llr_prior)
    else
        if length(initial_marginals) != num_bits
            error("initial_marginals has to have length equal to the number of columns of H.")
        end
        marginals_prev = Float64.(collect(initial_marginals))
    end

    gamma_vec =
        if gamma isa Real
            fill(Float64(gamma), num_bits)
        else
            if length(gamma) != num_bits
                error("gamma has to be a scalar or a vector of length equal to the number of columns of H.")
            end
            Float64.(collect(gamma))
        end

    check_to_bits, bit_to_checks, check_bit_pos, bit_check_pos = _relay_build_graph(H)

    var_to_check = [fill(llr_prior[j], length(bit_to_checks[j])) for j in 1:num_bits]
    check_to_var = [zeros(Float64, length(check_to_bits[i])) for i in 1:num_checks]

    bias = copy(llr_prior)
    marginals = copy(marginals_prev)
    hard_error = zeros(Int64, num_bits)

    _bp_sign(x::Real) = x < 0 ? -1.0 : 1.0

    for iter in 1:max_iter
        bias = (1 .- gamma_vec) .* llr_prior .+ gamma_vec .* marginals_prev

        for check_idx in 1:num_checks
            neighbors = check_to_bits[check_idx]
            deg = length(neighbors)
            if deg == 0
                continue
            end

            signs = Vector{Float64}(undef, deg)
            absvals = Vector{Float64}(undef, deg)
            total_sign = 1.0
            min1 = Inf
            min2 = Inf
            min1_idx = 0

            for local_idx in 1:deg
                bit_idx = neighbors[local_idx]
                msg = var_to_check[bit_idx][bit_check_pos[bit_idx][check_idx]]
                signs[local_idx] = _bp_sign(msg)
                absvals[local_idx] = abs(msg)
                total_sign *= signs[local_idx]

                if absvals[local_idx] < min1
                    min2 = min1
                    min1 = absvals[local_idx]
                    min1_idx = local_idx
                elseif absvals[local_idx] < min2
                    min2 = absvals[local_idx]
                end
            end

            syndrome_sign = s[check_idx] == 0 ? 1.0 : -1.0
            for local_idx in 1:deg
                excluded_sign = total_sign * signs[local_idx]
                min_without_j = local_idx == min1_idx ? min2 : min1
                if isinf(min_without_j)
                    min_without_j = 0.0
                end
                check_to_var[check_idx][local_idx] = syndrome_sign * excluded_sign * min_without_j
            end
        end

        for bit_idx in 1:num_bits
            neighbors = bit_to_checks[bit_idx]
            deg = length(neighbors)
            if deg == 0
                marginals[bit_idx] = bias[bit_idx]
                hard_error[bit_idx] = bias[bit_idx] < 0 ? 1 : 0
                continue
            end

            total = bias[bit_idx]
            for check_idx in neighbors
                total += check_to_var[check_idx][check_bit_pos[check_idx][bit_idx]]
            end
            marginals[bit_idx] = total
            hard_error[bit_idx] = total < 0 ? 1 : 0

            for local_idx in 1:deg
                check_idx = neighbors[local_idx]
                msg = bias[bit_idx]
                for other_check_idx in neighbors
                    if other_check_idx != check_idx
                        msg += check_to_var[other_check_idx][check_bit_pos[other_check_idx][bit_idx]]
                    end
                end
                var_to_check[bit_idx][local_idx] = msg
            end
        end

        syndrome_residual = mod.(H * hard_error .+ s, 2)
        if all(x -> x == 0, syndrome_residual)
            return (
                converged=true,
                error=copy(hard_error),
                marginals=copy(marginals),
                bias=copy(bias),
                iterations=iter,
                syndrome_weight=0,
                weight=_relay_weight(hard_error, llr_prior),
                gamma=copy(gamma_vec),
            )
        end

        marginals_prev = copy(marginals)
    end

    syndrome_residual = mod.(H * hard_error .+ s, 2)
    return (
        converged=false,
        error=copy(hard_error),
        marginals=copy(marginals),
        bias=copy(bias),
        iterations=max_iter,
        syndrome_weight=count(x -> x != 0, syndrome_residual),
        weight=_relay_weight(hard_error, llr_prior),
        gamma=copy(gamma_vec),
    )
end

"""
    relay_bp_decode(H::AbstractMatrix{<:Integer}, s::AbstractVector{<:Integer}, p;
                    num_solutions::Int=1,
                    max_legs::Int=301,
                    leg_max_iter::Int=60,
                    first_leg_max_iter::Int=80,
                    first_leg_gamma::Real=0.35,
                    gamma_center::Real=0.3655,
                    gamma_width::Real=1.239,
                    gamma_schedule=nothing,
                    rng::AbstractRNG=Random.default_rng())

Decode the binary syndrome equation `H * e = s (mod 2)` using the Relay-BP
heuristic introduced by Müller et al. (IBM Quantum).

Relay-BP chains together multiple DMem-BP legs. The first leg starts from the
channel log-likelihood ratios. Each subsequent leg is initialized with the
previous leg's final marginals, and uses a new set of memory strengths.
Whenever a leg finds a syndrome-matching solution, that solution is recorded.
The decoder stops once either `max_legs` legs have been executed or
`num_solutions` distinct solutions have been found, and returns the
lowest-weight converged solution.

The default memory-strength parameters are the rotated-surface-code XZ-decoding
settings reported by Müller et al.: first-leg memory strength `0.35`, followed
by i.i.d. sampling from the interval `[-0.254, 0.985]` on later legs.

Keyword arguments:
- `num_solutions`: stop after this many distinct converged solutions are found
- `max_legs`: maximum number of relay legs
- `leg_max_iter`: iteration limit for legs `2, 3, ...`
- `first_leg_max_iter`: iteration limit for the first leg
- `first_leg_gamma`: uniform memory strength used on the first leg when
  `gamma_schedule` is not supplied
- `gamma_center`, `gamma_width`: define the later-leg sampling interval
  `[gamma_center - gamma_width / 2, gamma_center + gamma_width / 2]`
- `gamma_schedule`: optional explicit schedule of gamma vectors; this can be
  either a matrix whose rows are gamma vectors or a vector of gamma vectors
- `rng`: random-number generator used for sampling gamma vectors

The return value is a named tuple with fields
- `error`: the returned correction (best converged correction if one exists,
  otherwise the best non-converged hard decision by residual-syndrome weight)
- `converged`: whether any leg converged
- `iterations`: the total number of BP iterations across all legs
- `num_legs_run`: the number of legs executed
- `solutions_found`: the number of distinct converged solutions found
- `weight`: the weighted cost of the returned correction
- `legs`: the per-leg decoder outputs
- `best_leg`: the index of the leg that produced the returned correction
"""
function relay_bp_decode(
    H::AbstractMatrix{<:Integer},
    s::AbstractVector{<:Integer},
    p;
    num_solutions::Int=1,
    max_legs::Int=301,
    leg_max_iter::Int=60,
    first_leg_max_iter::Int=80,
    first_leg_gamma::Real=0.35,
    gamma_center::Real=0.3655,
    gamma_width::Real=1.239,
    gamma_schedule=nothing,
    rng::AbstractRNG=Random.default_rng(),
)
    H = mod.(Int64.(H), 2)
    s = mod.(Int64.(s), 2)

    num_checks, num_bits = size(H)
    if length(s) != num_checks
        error("The syndrome length has to equal the number of rows of H.")
    end
    if num_solutions < 1
        error("num_solutions has to be a positive integer.")
    end
    if max_legs < 1
        error("max_legs has to be a positive integer.")
    end
    if leg_max_iter < 1 || first_leg_max_iter < 1
        error("Iteration limits have to be positive integers.")
    end

    llr_prior = _relay_prior_llr_vector(p, num_bits)
    initial_marginals = copy(llr_prior)

    legs = NamedTuple[]
    total_iterations = 0

    best_solution = nothing
    best_solution_leg = 0
    best_fallback = nothing
    best_fallback_leg = 0

    found_solution_keys = Set{BitVector}()

    for leg_idx in 1:max_legs
        gamma_vec = _relay_sample_gamma_vector(
            rng,
            num_bits,
            leg_idx;
            gamma_schedule=gamma_schedule,
            first_leg_gamma=first_leg_gamma,
            gamma_center=gamma_center,
            gamma_width=gamma_width,
        )

        max_iter = leg_idx == 1 ? first_leg_max_iter : leg_max_iter

        leg = dmem_bp_decode(
            H,
            s,
            p;
            max_iter=max_iter,
            initial_marginals=initial_marginals,
            gamma=gamma_vec,
        )

        push!(legs, merge((leg=leg_idx,), leg))
        total_iterations += leg.iterations
        initial_marginals = leg.marginals

        if isnothing(best_fallback)
            best_fallback = leg
            best_fallback_leg = leg_idx
        else
            better_residual = leg.syndrome_weight < best_fallback.syndrome_weight
            same_residual_better_weight =
                leg.syndrome_weight == best_fallback.syndrome_weight && leg.weight < best_fallback.weight
            if better_residual || same_residual_better_weight
                best_fallback = leg
                best_fallback_leg = leg_idx
            end
        end

        if leg.converged
            key = BitVector(leg.error .!= 0)
            if !(key in found_solution_keys)
                push!(found_solution_keys, key)
                if isnothing(best_solution) || leg.weight < best_solution.weight
                    best_solution = leg
                    best_solution_leg = leg_idx
                end
                if length(found_solution_keys) >= num_solutions
                    break
                end
            end
        end
    end

    if !isnothing(best_solution)
        return (
            error=best_solution.error,
            converged=true,
            iterations=total_iterations,
            num_legs_run=length(legs),
            solutions_found=length(found_solution_keys),
            weight=best_solution.weight,
            legs=legs,
            best_leg=best_solution_leg,
        )
    else
        return (
            error=best_fallback.error,
            converged=false,
            iterations=total_iterations,
            num_legs_run=length(legs),
            solutions_found=0,
            weight=best_fallback.weight,
            legs=legs,
            best_leg=best_fallback_leg,
        )
    end
end

"""
    css_relay_bp_decode(HX::AbstractMatrix{<:Integer}, HZ::AbstractMatrix{<:Integer},
                        sx::AbstractVector{<:Integer}, sz::AbstractVector{<:Integer}, p;
                        kwargs...)

Decode a CSS code under uncorrelated code-capacity `X/Z` noise using two
independent calls to `relay_bp_decode`.

The convention is:
- `sx = HZ * x (mod 2)` is the syndrome induced by `X` errors
- `sz = HX * z (mod 2)` is the syndrome induced by `Z` errors

All keyword arguments are passed through to `relay_bp_decode` in each sector.
The return value is a named tuple with fields
- `x`: the estimated `X`-component correction
- `z`: the estimated `Z`-component correction
- `x_result`: the full decoder output for the `X` sector
- `z_result`: the full decoder output for the `Z` sector
"""
function css_relay_bp_decode(
    HX::AbstractMatrix{<:Integer},
    HZ::AbstractMatrix{<:Integer},
    sx::AbstractVector{<:Integer},
    sz::AbstractVector{<:Integer},
    p;
    kwargs...
)
    x_result = relay_bp_decode(HZ, sx, p; kwargs...)
    z_result = relay_bp_decode(HX, sz, p; kwargs...)
    return (x=x_result.error, z=z_result.error, x_result=x_result, z_result=z_result)
end
