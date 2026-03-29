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

function _bp_prepare_initial_marginals(initial_marginals, llr_prior::Vector{Float64}, num_bits::Int)
    if isnothing(initial_marginals)
        return copy(llr_prior)
    end
    if length(initial_marginals) != num_bits
        error("initial_marginals has to have length equal to the number of columns of H.")
    end
    return Float64.(collect(initial_marginals))
end

function _bp_precompute_message_positions(
    check_to_bits::Vector{Vector{Int64}},
    bit_to_checks::Vector{Vector{Int64}},
    check_bit_pos,
    bit_check_pos,
)
    check_neighbor_pos_in_bit = Vector{Vector{Int64}}(undef, length(check_to_bits))
    for check_idx in eachindex(check_to_bits)
        neighbors = check_to_bits[check_idx]
        positions = Vector{Int64}(undef, length(neighbors))
        @inbounds for local_idx in 1:length(neighbors)
            bit_idx = neighbors[local_idx]
            positions[local_idx] = bit_check_pos[bit_idx][check_idx]
        end
        check_neighbor_pos_in_bit[check_idx] = positions
    end

    bit_neighbor_pos_in_check = Vector{Vector{Int64}}(undef, length(bit_to_checks))
    for bit_idx in eachindex(bit_to_checks)
        neighbors = bit_to_checks[bit_idx]
        positions = Vector{Int64}(undef, length(neighbors))
        @inbounds for local_idx in 1:length(neighbors)
            check_idx = neighbors[local_idx]
            positions[local_idx] = check_bit_pos[check_idx][bit_idx]
        end
        bit_neighbor_pos_in_check[bit_idx] = positions
    end

    return check_neighbor_pos_in_bit, bit_neighbor_pos_in_check
end

function _bp_fill_bias!(
    bias::Vector{Float64},
    llr_prior::Vector{Float64},
    γ_vec::Vector{Float64},
    marginals_prev::Vector{Float64},
    bit_to_check_update_rule::Symbol,
)
    if bit_to_check_update_rule == :memoryless
        return nothing
    end

    @inbounds for j in eachindex(bias)
        bias[j] = (1 - γ_vec[j]) * llr_prior[j] + γ_vec[j] * marginals_prev[j]
    end
    return nothing
end

function _bp_weight(hard_error::Vector{Int64}, llr_prior::Vector{Float64})
    total = 0.0
    @inbounds for j in eachindex(hard_error)
        total += hard_error[j] * llr_prior[j]
    end
    return total
end

function _bp_syndrome_weight(
    check_to_bits::Vector{Vector{Int64}},
    hard_error::Vector{Int64},
    s::AbstractVector{<:Integer},
)
    residual_weight = 0
    @inbounds for check_idx in eachindex(check_to_bits)
        parity = s[check_idx]
        for bit_idx in check_to_bits[check_idx]
            parity = xor(parity, hard_error[bit_idx])
        end
        residual_weight += parity
    end
    return residual_weight
end

function _bp_check_update_sum_product_from_graph!(
    outgoing_messages::Vector{Float64},
    neighbors::Vector{Int64},
    neighbor_positions_in_bit_msgs::Vector{Int64},
    bit_to_check,
    syndrome_sign::Float64,
    tanh_half_workspace::Vector{Float64};
    α::Real=1.0,
)
    deg = length(neighbors)
    if deg == 0
        return nothing
    end

    if deg == 1
        sat = 2 * atanh(1 - 1e-15)
        outgoing_messages[1] = Float64(α) * syndrome_sign * sat
        return nothing
    end

    zero_count = 0
    zero_idx = 0
    product_nonzero = 1.0

    @inbounds for local_idx in 1:deg
        bit_idx = neighbors[local_idx]
        msg = bit_to_check[bit_idx][neighbor_positions_in_bit_msgs[local_idx]]
        th = tanh(msg / 2)
        tanh_half_workspace[local_idx] = th
        if iszero(th)
            zero_count += 1
            zero_idx = local_idx
        else
            product_nonzero *= th
        end
    end

    scale = Float64(α) * syndrome_sign * 2.0

    if zero_count > 1
        @inbounds for local_idx in 1:deg
            outgoing_messages[local_idx] = 0.0
        end
    elseif zero_count == 1
        @inbounds for local_idx in 1:deg
            prod_excluding_j = local_idx == zero_idx ? product_nonzero : 0.0
            clipped = clamp(prod_excluding_j, -1 + 1e-15, 1 - 1e-15)
            outgoing_messages[local_idx] = scale * atanh(clipped)
        end
    else
        @inbounds for local_idx in 1:deg
            prod_excluding_j = product_nonzero / tanh_half_workspace[local_idx]
            clipped = clamp(prod_excluding_j, -1 + 1e-15, 1 - 1e-15)
            outgoing_messages[local_idx] = scale * atanh(clipped)
        end
    end

    return nothing
end

function _bp_check_update_min_sum_from_graph!(
    outgoing_messages::Vector{Float64},
    neighbors::Vector{Int64},
    neighbor_positions_in_bit_msgs::Vector{Int64},
    bit_to_check,
    syndrome_sign::Float64;
    α::Real=1.0,
)
    deg = length(neighbors)
    if deg == 0
        return nothing
    end

    if deg == 1
        sat = 2 * atanh(1 - 1e-15)
        outgoing_messages[1] = Float64(α) * syndrome_sign * sat
        return nothing
    end

    total_sign = 1.0
    min1 = Inf
    min2 = Inf
    min1_idx = 0

    @inbounds for local_idx in 1:deg
        bit_idx = neighbors[local_idx]
        msg = bit_to_check[bit_idx][neighbor_positions_in_bit_msgs[local_idx]]
        sign = msg < 0 ? -1.0 : 1.0
        absmsg = abs(msg)
        total_sign *= sign

        if absmsg < min1
            min2 = min1
            min1 = absmsg
            min1_idx = local_idx
        elseif absmsg < min2
            min2 = absmsg
        end
    end

    scale = syndrome_sign * Float64(α)
    @inbounds for local_idx in 1:deg
        bit_idx = neighbors[local_idx]
        msg = bit_to_check[bit_idx][neighbor_positions_in_bit_msgs[local_idx]]
        sign = msg < 0 ? -1.0 : 1.0
        excluded_sign = total_sign * sign
        min_without_j = local_idx == min1_idx ? min2 : min1
        outgoing_messages[local_idx] = scale * excluded_sign * min_without_j
    end

    return nothing
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
    initial_beliefs = _bp_prepare_initial_marginals(initial_marginals, llr_prior, num_bits)

    check_to_bits, bit_to_checks, check_bit_pos, bit_check_pos = tanner_graph(H)
    check_neighbor_pos_in_bit, bit_neighbor_pos_in_check = _bp_precompute_message_positions(
        check_to_bits,
        bit_to_checks,
        check_bit_pos,
        bit_check_pos,
    )

    bit_to_check = Vector{Vector{Float64}}(undef, num_bits)
    @inbounds for bit_idx in 1:num_bits
        bit_to_check[bit_idx] = fill(initial_beliefs[bit_idx], length(bit_to_checks[bit_idx]))
    end
    check_to_bit = [zeros(Float64, length(check_to_bits[check_idx])) for check_idx in 1:num_checks]

    bias = copy(llr_prior)
    marginals_prev = copy(initial_beliefs)
    marginals = similar(initial_beliefs)
    hard_error = zeros(Int64, num_bits)

    is_memoryless = bit_to_check_update_rule == :memoryless
    use_sum_product = check_to_bit_update_rule == :sum_product
    max_check_degree = num_checks == 0 ? 0 : maximum(length(check_to_bits[check_idx]) for check_idx in 1:num_checks)
    tanh_half_workspace = use_sum_product ? zeros(Float64, max_check_degree) : Float64[]

    for iter in 1:max_iter
        if !is_memoryless
            _bp_fill_bias!(bias, llr_prior, γ_vec, marginals_prev, bit_to_check_update_rule)
        end

        α = if check_to_bit_update_rule == :min_sum
            min_sum_scaling == :roffe ? 1 - 2.0^(-iter) : 1.0
        else
            1.0
        end

        if use_sum_product
            @inbounds for check_idx in 1:num_checks
                syndrome_sign = s[check_idx] == 0 ? 1.0 : -1.0
                _bp_check_update_sum_product_from_graph!(
                    check_to_bit[check_idx],
                    check_to_bits[check_idx],
                    check_neighbor_pos_in_bit[check_idx],
                    bit_to_check,
                    syndrome_sign,
                    tanh_half_workspace;
                    α=α,
                )
            end
        else
            @inbounds for check_idx in 1:num_checks
                syndrome_sign = s[check_idx] == 0 ? 1.0 : -1.0
                _bp_check_update_min_sum_from_graph!(
                    check_to_bit[check_idx],
                    check_to_bits[check_idx],
                    check_neighbor_pos_in_bit[check_idx],
                    bit_to_check,
                    syndrome_sign;
                    α=α,
                )
            end
        end

        @inbounds for bit_idx in 1:num_bits
            neighbors = bit_to_checks[bit_idx]
            incoming_positions = bit_neighbor_pos_in_check[bit_idx]
            deg = length(neighbors)

            total = bias[bit_idx]
            for local_idx in 1:deg
                check_idx = neighbors[local_idx]
                total += check_to_bit[check_idx][incoming_positions[local_idx]]
            end

            marginals[bit_idx] = total
            hard_error[bit_idx] = total < 0 ? 1 : 0

            for local_idx in 1:deg
                check_idx = neighbors[local_idx]
                bit_to_check[bit_idx][local_idx] = total - check_to_bit[check_idx][incoming_positions[local_idx]]
            end
        end

        syndrome_weight = _bp_syndrome_weight(check_to_bits, hard_error, s)
        if syndrome_weight == 0
            return (
                converged=true,
                error=copy(hard_error),
                marginals=copy(marginals),
                llr=copy(marginals),
                bias=copy(bias),
                iterations=iter,
                syndrome_weight=0,
                weight=_bp_weight(hard_error, llr_prior),
                check_to_bit_update_rule=check_to_bit_update_rule,
                bit_to_check_update_rule=bit_to_check_update_rule,
                gamma=copy(γ_vec),
            )
        end

        if !is_memoryless
            marginals_prev, marginals = marginals, marginals_prev
        end
    end

    final_marginals = is_memoryless ? marginals : marginals_prev
    syndrome_weight = _bp_syndrome_weight(check_to_bits, hard_error, s)
    return (
        converged=false,
        error=copy(hard_error),
        marginals=copy(final_marginals),
        llr=copy(final_marginals),
        bias=copy(bias),
        iterations=max_iter,
        syndrome_weight=syndrome_weight,
        weight=_bp_weight(hard_error, llr_prior),
        check_to_bit_update_rule=check_to_bit_update_rule,
        bit_to_check_update_rule=bit_to_check_update_rule,
        gamma=copy(γ_vec),
    )
end
