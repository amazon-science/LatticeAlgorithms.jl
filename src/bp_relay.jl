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

"""
    _relay_sample_gamma_vector(rng::AbstractRNG, num_bits::Int, leg_idx::Int;
                               gamma_schedule=nothing,
                               first_leg_gamma::Real=0.35,
                               gamma_center::Real=0.3655,
                               gamma_width::Real=1.239)

Return the length-`num_bits` memory-strength vector used by one relay leg.

If `gamma_schedule` is supplied, it overrides the default sampling rule.
Supported formats are

- a matrix whose `leg_idx`-th row is the gamma vector for that leg;
- a vector whose `leg_idx`-th entry is either
  - a scalar gamma value, or
  - a length-`num_bits` gamma vector.

If `gamma_schedule` is not supplied, the default rule is

- first leg: a uniform vector with value `first_leg_gamma`;
- later legs: i.i.d. sampling from
  `[gamma_center - gamma_width / 2, gamma_center + gamma_width / 2]`.
"""
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

            γ_leg = gamma_schedule[leg_idx]

            if γ_leg isa Real
                γf = Float64(γ_leg)
                isfinite(γf) || error("Scalar entries of gamma_schedule must be finite.")
                return fill(γf, num_bits)
            else
                if length(γ_leg) != num_bits
                    error("Each gamma vector in gamma_schedule has to have length num_bits.")
                end
                γ_vec = Float64.(collect(γ_leg))
                all(isfinite, γ_vec) || error("All entries of gamma_schedule must be finite.")
                return γ_vec
            end
        else
            error("gamma_schedule has to be a matrix or a vector whose entries are scalars or vectors.")
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

"""
    relay_bp_decode(H::AbstractMatrix{<:Integer},
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
                    check_to_bit_update_rule::Symbol=:min_sum,
                    min_sum_scaling::Symbol=:none,
                    rng::AbstractRNG=Random.default_rng())

Decode the binary syndrome equation `H * e = s (mod 2)` using the Relay-BP
heuristic.

Each relay leg is a call to

```julia
bp_decode(...;
    check_to_bit_update_rule=check_to_bit_update_rule,
    bit_to_check_update_rule=:dmem,
    γ=γ_leg,
    initial_marginals=previous_leg_marginals,
    min_sum_scaling=min_sum_scaling)
```

So this file only implements the relay scheduler; the underlying BP dynamics are
delegated to `bp_decode` in `bp.jl`.

The canonical Relay-BP setting is
- `check_to_bit_update_rule = :min_sum`
- `min_sum_scaling = :none`

The first leg starts from the channel prior. Each subsequent leg is initialized
with the previous leg's final marginals and uses a new gamma vector.

Whenever a leg finds a syndrome-consistent solution, that solution is recorded.
The decoder stops once either
- `max_legs` legs have been executed, or
- `num_solutions` distinct converged solutions have been found.

Among converged solutions, the decoder returns the one of minimum weighted cost.
If no leg converges, it returns the best fallback hard decision, ranked first by
residual syndrome weight and then by weighted cost.

Keyword arguments
-----------------
- `num_solutions`: stop after this many distinct converged solutions are found.
- `max_legs`: maximum number of relay legs.
- `leg_max_iter`: iteration limit for legs `2, 3, ...`.
- `first_leg_max_iter`: iteration limit for the first leg.
- `first_leg_gamma`: uniform memory strength used on the first leg when
  `gamma_schedule` is not supplied.
- `gamma_center`, `gamma_width`: define the later-leg sampling interval
  `[gamma_center - gamma_width / 2, gamma_center + gamma_width / 2]`.
- `gamma_schedule`: optional explicit schedule of gamma values/vectors.
- `check_to_bit_update_rule`: passed through to `bp_decode`.
- `min_sum_scaling`: passed through to `bp_decode`.
- `rng`: random-number generator used for sampling gamma vectors.

Return value
------------
A named tuple with fields
- `error`
- `converged`
- `iterations`
- `num_legs_run`
- `solutions_found`
- `weight`
- `legs`
- `best_leg`
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
    check_to_bit_update_rule::Symbol=:min_sum,
    min_sum_scaling::Symbol=:none,
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

    initial_marginals = nothing
    legs = NamedTuple[]
    total_iterations = 0

    best_solution = nothing
    best_solution_leg = 0
    best_fallback = nothing
    best_fallback_leg = 0

    found_solution_keys = Set{BitVector}()

    for leg_idx in 1:max_legs
        γ_vec = _relay_sample_gamma_vector(
            rng,
            num_bits,
            leg_idx;
            gamma_schedule=gamma_schedule,
            first_leg_gamma=first_leg_gamma,
            gamma_center=gamma_center,
            gamma_width=gamma_width,
        )

        max_iter = leg_idx == 1 ? first_leg_max_iter : leg_max_iter

        leg = bp_decode(
            H,
            s,
            p;
            max_iter=max_iter,
            check_to_bit_update_rule=check_to_bit_update_rule,
            bit_to_check_update_rule=:dmem,
            γ=γ_vec,
            initial_marginals=initial_marginals,
            min_sum_scaling=min_sum_scaling,
        )

        push!(legs, merge((leg=leg_idx,), leg))
        total_iterations += leg.iterations

        # Relay initialization: the next leg starts from this leg's final marginals.
        initial_marginals = leg.marginals

        # Track the best non-converged fallback by residual syndrome weight, then by cost.
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

        # Record distinct converged solutions.
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
    css_relay_bp_decode(HX::AbstractMatrix{<:Integer},
                        HZ::AbstractMatrix{<:Integer},
                        sx::AbstractVector{<:Integer},
                        sz::AbstractVector{<:Integer},
                        px,
                        pz=px;
                        kwargs...)

Decode a CSS code under independent code-capacity `X/Z` noise using two
independent calls to `relay_bp_decode`.

The convention is
- `sx = HZ * x (mod 2)` for the `X`-error sector,
- `sz = HX * z (mod 2)` for the `Z`-error sector.

All keyword arguments are passed through to `relay_bp_decode` in each sector.

Return value
------------
A named tuple with fields
- `x`
- `z`
- `x_result`
- `z_result`
"""
function css_relay_bp_decode(
    HX::AbstractMatrix{<:Integer},
    HZ::AbstractMatrix{<:Integer},
    sx::AbstractVector{<:Integer},
    sz::AbstractVector{<:Integer},
    px,
    pz=px;
    kwargs...,
)
    x_result = relay_bp_decode(HZ, sx, px; kwargs...)
    z_result = relay_bp_decode(HX, sz, pz; kwargs...)

    return (
        x=x_result.error,
        z=z_result.error,
        x_result=x_result,
        z_result=z_result,
    )
end
