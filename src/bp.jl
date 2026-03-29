# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# Licensed under the Apache License, Version 2.0 (the "License").
# You may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

using SparseArrays

"""
    _bp_prior_llr_vector(p::Real, num_bits::Int)
    _bp_prior_llr_vector(p::AbstractVector{<:Real}, num_bits::Int)

Return the channel log-likelihood-ratio vector

```math
Lambda_j^{(0)} = log(frac{1-p_j}{p_j})
```

for a binary symmetric channel. The argument `p` can be either a scalar bit-flip
probability or a vector of per-bit probabilities.
"""
function _bp_prior_llr_vector(p::Real, num_bits::Int)
    if !(0 < p < 0.5)
        error("The bit-flip probability p has to satisfy 0 < p < 0.5.")
    end
    return fill(log((1 - p) / p), num_bits)
end

function _bp_prior_llr_vector(p::AbstractVector{<:Real}, num_bits::Int)
    if length(p) != num_bits
        error("The probability vector length has to equal the number of columns of H.")
    end
    if any(x -> !(0 < x < 0.5), p)
        error("All bit-flip probabilities have to satisfy 0 < p_j < 0.5.")
    end
    return log.((1 .- p) ./ p)
end

"""
    _bp_validate_update_rules(check_to_bit_update_rule::Symbol,
                              bit_to_check_update_rule::Symbol)

Validate the two BP update-rule selectors.
"""
function _bp_validate_update_rules(
    check_to_bit_update_rule::Symbol,
    bit_to_check_update_rule::Symbol,
)
    if check_to_bit_update_rule ∉ (:sum_product, :min_sum)
        error(
            "Unsupported check_to_bit_update_rule=$(check_to_bit_update_rule). " *
            "Supported values are :sum_product and :min_sum.",
        )
    end
    if bit_to_check_update_rule ∉ (:memoryless, :mem, :dmem)
        error(
            "Unsupported bit_to_check_update_rule=$(bit_to_check_update_rule). " *
            "Supported values are :memoryless, :mem, and :dmem.",
        )
    end
    return nothing
end

"""
    _bp_gamma_vector(γ, bit_to_check_update_rule::Symbol, num_bits::Int)

Normalize the memory parameter to a length-`num_bits` vector.
"""
function _bp_gamma_vector(γ, bit_to_check_update_rule::Symbol, num_bits::Int)
    if bit_to_check_update_rule == :memoryless
        if isnothing(γ)
            return zeros(Float64, num_bits)
        elseif γ isa Real
            if γ != 0
                error("For bit_to_check_update_rule=:memoryless, γ must be omitted or set to 0.")
            end
            return zeros(Float64, num_bits)
        else
            if length(γ) != num_bits
                error(
                    "For bit_to_check_update_rule=:memoryless, a vector γ must have length num_bits " *
                    "and contain only zeros.",
                )
            end
            γ_vec = Float64.(collect(γ))
            if any(x -> !iszero(x), γ_vec)
                error("For bit_to_check_update_rule=:memoryless, γ must be identically zero.")
            end
            return γ_vec
        end
    elseif bit_to_check_update_rule == :mem
        if isnothing(γ)
            error("For bit_to_check_update_rule=:mem, γ must be supplied as a scalar.")
        elseif γ isa Real
            γf = Float64(γ)
            if !(0.0 <= γf <= 1.0)
                error("For bit_to_check_update_rule=:mem, γ must be a scalar in [0, 1].")
            end
            return fill(γf, num_bits)
        else
            error("For bit_to_check_update_rule=:mem, γ must be a scalar.")
        end
    else # :dmem
        if isnothing(γ)
            error("For bit_to_check_update_rule=:dmem, γ must be supplied as a scalar or a vector.")
        elseif γ isa Real
            return fill(Float64(γ), num_bits)
        else
            if length(γ) != num_bits
                error(
                    "For bit_to_check_update_rule=:dmem, γ must be a scalar or a vector of length num_bits.",
                )
            end
            return Float64.(collect(γ))
        end
    end
end

"""
    _bp_check_update_sum_product!(outgoing_messages, incoming_messages, syndrome_sign; α=1.0)

Apply the essential LLR-domain sum-product local rule for one check node.
The caller is responsible for gathering the incoming messages around the check.
The optional scaling factor `α` multiplies the final outgoing message.
"""
function _bp_check_update_sum_product!(
    outgoing_messages::Vector{Float64},
    incoming_messages::Vector{Float64},
    syndrome_sign::Float64;
    α::Real=1.0,
)
    deg = length(incoming_messages)
    tanh_half = Vector{Float64}(undef, deg)

    for local_idx in 1:deg
        tanh_half[local_idx] = tanh(incoming_messages[local_idx] / 2)
    end

    total_prod = prod(tanh_half)

    for local_idx in 1:deg
        x = tanh_half[local_idx]
        prod_excluding_j = iszero(x) ? prod(tanh_half[k] for k in 1:deg if k != local_idx) : total_prod / x

        # Guard against roundoff slightly outside (-1, 1).
        clipped = clamp(prod_excluding_j, -1 + 1e-15, 1 - 1e-15)
        outgoing_messages[local_idx] = Float64(α) * syndrome_sign * 2 * atanh(clipped)
    end
end

"""
    _bp_check_update_min_sum!(outgoing_messages, incoming_messages, syndrome_sign; α=1.0)

Apply the essential min-sum local rule for one check node. The caller is
responsible for gathering the incoming messages around the check.
"""
function _bp_check_update_min_sum!(
    outgoing_messages::Vector{Float64},
    incoming_messages::Vector{Float64},
    syndrome_sign::Float64;
    α::Real=1.0,
)
    deg = length(incoming_messages)

    # Degree-1 check: the check fixes the bit directly, so send a saturated LLR.
    if deg == 1
        sat = 2 * atanh(1 - 1e-15)
        outgoing_messages[1] = Float64(α) * syndrome_sign * sat
        return nothing
    end

    signs = Vector{Float64}(undef, deg)
    absvals = Vector{Float64}(undef, deg)
    total_sign = 1.0
    min1 = Inf
    min2 = Inf
    min1_idx = 0

    for local_idx in 1:deg
        msg = incoming_messages[local_idx]
        signs[local_idx] = msg < 0 ? -1.0 : 1.0
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

    for local_idx in 1:deg
        excluded_sign = total_sign * signs[local_idx]
        min_without_j = local_idx == min1_idx ? min2 : min1
        outgoing_messages[local_idx] = syndrome_sign * Float64(α) * excluded_sign * min_without_j
    end
end


"""
    bp_decode(H::AbstractMatrix{<:Integer}, s::AbstractVector{<:Integer}, p;
              max_iter::Int=size(H, 2),
              check_to_bit_update_rule::Symbol=:sum_product,
              bit_to_check_update_rule::Symbol=:memoryless,
              γ=nothing,
              initial_marginals=nothing,
              min_sum_scaling::Symbol=:none)

Decode the binary syndrome equation `H * e = s (mod 2)` with a configurable
belief-propagation engine.

The implementation exposes two orthogonal update-rule choices.

- `check_to_bit_update_rule`:
  - `:sum_product` for the exact LLR-domain tanh/atanh check update,
  - `:min_sum` for the min-sum approximation.
- `bit_to_check_update_rule`:
  - `:memoryless` for the standard BP variable-node update,
  - `:mem` for uniform-memory BP,
  - `:dmem` for disordered-memory BP.

This means the common named decoders correspond to the following keyword pairs.

- Standard BP (sum-product / memoryless BP):
  `check_to_bit_update_rule=:sum_product`,
  `bit_to_check_update_rule=:memoryless`.
- Min-sum BP:
  `check_to_bit_update_rule=:min_sum`,
  `bit_to_check_update_rule=:memoryless`.
- Mem-BP:
  canonically `check_to_bit_update_rule=:sum_product`,
  `bit_to_check_update_rule=:mem`, with scalar `γ` between 0 and 1.
- DMem-BP:
  canonically `check_to_bit_update_rule=:sum_product` or `:min_sum`,
  `bit_to_check_update_rule=:dmem`, with node-dependent `γ`.
  In Relay-BP applications, one often combines `:dmem` with `:min_sum`.

The channel parameter `p` can be either a scalar bit-flip probability or a
vector of per-bit probabilities.

Keyword arguments
-----------------
- `max_iter`: maximum number of BP iterations.
- `γ`: memory parameter. Use a scalar for `:mem`, and either a scalar or a
  length-`n` vector for `:dmem`.
- `initial_marginals`: initial beliefs `B^(0)`. If omitted, the channel LLRs are
  used.
- `min_sum_scaling`:
  - `:none` for unscaled min-sum,
  - `:roffe` for the schedule `α_t = 1 - 2^(-t)` used in the current
    min-sum BP implementation.

Return value
------------
A named tuple with fields
- `converged`
- `error`
- `marginals`
- `llr` (an alias of `marginals`, for compatibility with BP+OSD code)
- `bias`
- `iterations`
- `syndrome_weight`
- `weight`
- `check_to_bit_update_rule`
- `bit_to_check_update_rule`
- `gamma`
"""
function bp_decode(
    H::AbstractMatrix{<:Integer},
    s::AbstractVector{<:Integer},
    p;
    max_iter::Int=size(H, 2),
    check_to_bit_update_rule::Symbol=:sum_product,
    bit_to_check_update_rule::Symbol=:memoryless,
    γ=nothing,
    initial_marginals=nothing,
    min_sum_scaling::Symbol=:none,
)
    _bp_validate_update_rules(check_to_bit_update_rule, bit_to_check_update_rule)

    if min_sum_scaling ∉ (:none, :roffe)
        error("Unsupported min_sum_scaling=$(min_sum_scaling). Supported values are :none and :roffe.")
    end

    H = sparse(mod.(Int64.(H), 2))
    s = mod.(Int64.(s), 2)
    num_checks, num_bits = size(H)

    if length(s) != num_checks
        error("The syndrome length has to equal the number of rows of H.")
    end
    if max_iter < 1
        error("max_iter has to be a positive integer.")
    end

    llr_prior = Float64.(_bp_prior_llr_vector(p, num_bits))
    γ_vec = _bp_gamma_vector(γ, bit_to_check_update_rule, num_bits)

    marginals_prev = if isnothing(initial_marginals)
        copy(llr_prior)
    else
        if length(initial_marginals) != num_bits
            error("initial_marginals has to have length equal to the number of columns of H.")
        end
        Float64.(collect(initial_marginals))
    end

    check_to_bits, bit_to_checks, check_bit_pos, bit_check_pos = tanner_graph(H)

    bit_to_check = [fill(marginals_prev[j], length(bit_to_checks[j])) for j in 1:num_bits]
    check_to_bit = [zeros(Float64, length(check_to_bits[i])) for i in 1:num_checks]

    marginals = copy(marginals_prev)
    hard_error = zeros(Int64, num_bits)
    for iter in 1:max_iter
        if bit_to_check_update_rule == :memoryless
            bias = copy(llr_prior)
        else
            bias = (1 .- γ_vec) .* llr_prior .+ γ_vec .* marginals_prev
        end

        α = if check_to_bit_update_rule == :min_sum
            min_sum_scaling == :roffe ? 1 - 2.0^(-iter) : 1.0
        else
            1.0
        end

        # Apply the check to bit message passing rule
        for check_idx in 1:num_checks
            neighbors = check_to_bits[check_idx]
            deg = length(neighbors)
            if deg == 0
                continue
            end

            # Gather all incoming bit-to-check messages attached to this check.
            incoming_messages = Vector{Float64}(undef, deg)
            for local_idx in 1:deg
                bit_idx = neighbors[local_idx]
                incoming_messages[local_idx] = bit_to_check[bit_idx][bit_check_pos[bit_idx][check_idx]]
            end

            syndrome_sign = s[check_idx] == 0 ? 1.0 : -1.0

            # Apply the selected local check rule to produce all outgoing check-to-bit messages.
            if check_to_bit_update_rule == :sum_product
                _bp_check_update_sum_product!(check_to_bit[check_idx], incoming_messages, syndrome_sign; α=α)
            else
                _bp_check_update_min_sum!(check_to_bit[check_idx], incoming_messages, syndrome_sign; α=α)
            end
        end

        # Apply the bit to check  message passing rule
        for bit_idx in 1:num_bits
            neighbors = bit_to_checks[bit_idx]
            deg = length(neighbors)

            if deg == 0
                marginals[bit_idx] = bias[bit_idx]
                hard_error[bit_idx] = bias[bit_idx] < 0 ? 1 : 0
                continue
            end

            # Combine the channel/memory bias with all incoming check-to-bit messages
            # to form the current marginal and hard decision on this bit.
            total = bias[bit_idx]
            for check_idx in neighbors
                total += check_to_bit[check_idx][check_bit_pos[check_idx][bit_idx]]
            end
            marginals[bit_idx] = total
            hard_error[bit_idx] = total < 0 ? 1 : 0

            # For each neighboring check, send the extrinsic bit-to-check message,
            # namely the bias plus all incoming check messages except the recipient one.
            for local_idx in 1:deg
                check_idx = neighbors[local_idx]
                msg = bias[bit_idx]
                for other_check_idx in neighbors
                    if other_check_idx != check_idx
                        msg += check_to_bit[other_check_idx][check_bit_pos[other_check_idx][bit_idx]]
                    end
                end
                bit_to_check[bit_idx][local_idx] = msg
            end
        end

        syndrome_residual = mod.(H * hard_error .+ s, 2)
        if all(iszero, syndrome_residual)
            return (
                converged=true,
                error=copy(hard_error),
                marginals=copy(marginals),
                llr=copy(marginals),
                bias=copy(bias),
                iterations=iter,
                syndrome_weight=0,
                weight=sum(Int64.(hard_error) .* llr_prior),
                check_to_bit_update_rule=check_to_bit_update_rule,
                bit_to_check_update_rule=bit_to_check_update_rule,
                gamma=copy(γ_vec),
            )
        end

        marginals_prev = copy(marginals)
    end

    syndrome_residual = mod.(H * hard_error .+ s, 2)
    return (
        converged=false,
        error=copy(hard_error),
        marginals=copy(marginals),
        llr=copy(marginals),
        bias=copy(bias),
        iterations=max_iter,
        syndrome_weight=count(x -> x != 0, syndrome_residual),
        weight=sum(Int64.(hard_error) .* llr_prior),
        check_to_bit_update_rule=check_to_bit_update_rule,
        bit_to_check_update_rule=bit_to_check_update_rule,
        gamma=copy(γ_vec),
    )
end
