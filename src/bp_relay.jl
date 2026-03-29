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
