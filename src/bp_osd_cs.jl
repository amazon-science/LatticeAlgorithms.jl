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

"""
    bp_min_sum_decode(H::AbstractMatrix{<:Integer}, s::AbstractVector{<:Integer}, p::Real;
                      max_iter::Int=size(H, 2))

Decode the binary syndrome equation `H * e = s (mod 2)` using the scaled
min-sum belief propagation algorithm described by
Roffe, White, Burton, and Campbell.

The channel is taken to be a binary symmetric channel with bit-flip rate `p`.

The return value is a named tuple with fields
- `converged`: whether BP found a candidate satisfying the syndrome equation,
- `error`: the BP hard decision,
- `llr`: the final soft-decision vector,
- `iterations`: the number of iterations performed.
"""
function bp_min_sum_decode(
    H::AbstractMatrix{<:Integer},
    s::AbstractVector{<:Integer},
    p::Real;
    max_iter::Int=size(H, 2),
)
    H = sparse(mod.(Int64.(H), 2))
    s = mod.(Int64.(s), 2)

    num_checks, num_bits = size(H)
    if length(s) != num_checks
        error("The syndrome length has to equal the number of rows of H.")
    end
    if !(0 < p < 0.5)
        error("The bit-flip probability p has to satisfy 0 < p < 0.5.")
    end
    if max_iter < 1
        error("max_iter has to be a positive integer.")
    end

    check_to_bits = [Int64[] for _ in 1:num_checks]
    bit_to_checks = [Int64[] for _ in 1:num_bits]

    # local position of bit_idx inside check_to_bits[check_idx]
    check_bit_pos = [Dict{Int64, Int64}() for _ in 1:num_checks]

    # local position of check_idx inside bit_to_checks[bit_idx]
    bit_check_pos = [Dict{Int64, Int64}() for _ in 1:num_bits]

    rows, cols, _ = findnz(H)
    for (check_idx, bit_idx) in zip(rows, cols)
        push!(check_to_bits[check_idx], bit_idx)
        check_bit_pos[check_idx][bit_idx] = length(check_to_bits[check_idx])

        push!(bit_to_checks[bit_idx], check_idx)
        bit_check_pos[bit_idx][check_idx] = length(bit_to_checks[bit_idx])
    end

    llr_channel = log((1 - p) / p)
    var_to_check = [fill(llr_channel, length(bit_to_checks[j])) for j in 1:num_bits]
    check_to_var = [zeros(Float64, length(check_to_bits[i])) for i in 1:num_checks]

    soft_llr = fill(llr_channel, num_bits)
    hard_error = zeros(Int64, num_bits)

    _bp_sign(x::Real) = x < 0 ? -1.0 : 1.0

    for iter in 1:max_iter
        α = 1 - 2.0^(-iter)

        # check -> variable updates
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
                check_to_var[check_idx][local_idx] = syndrome_sign * α * excluded_sign * min_without_j
            end
        end

        # variable -> check updates and hard decision
        for bit_idx in 1:num_bits
            neighbors = bit_to_checks[bit_idx]
            deg = length(neighbors)
            if deg == 0
                soft_llr[bit_idx] = llr_channel
                hard_error[bit_idx] = 0
                continue
            end

            total = llr_channel
            for check_idx in neighbors
                total += check_to_var[check_idx][check_bit_pos[check_idx][bit_idx]]
            end
            soft_llr[bit_idx] = total
            hard_error[bit_idx] = total < 0 ? 1 : 0

            for local_idx in 1:deg
                check_idx = neighbors[local_idx]
                msg = llr_channel
                for other_check_idx in neighbors
                    if other_check_idx != check_idx
                        msg += check_to_var[other_check_idx][check_bit_pos[other_check_idx][bit_idx]]
                    end
                end
                var_to_check[bit_idx][local_idx] = msg
            end
        end

        if mod.(H * hard_error, 2) == s
            return (converged=true, error=hard_error, llr=soft_llr, iterations=iter)
        end
    end

    return (converged=false, error=hard_error, llr=soft_llr, iterations=max_iter)
end

function _bp_osd_setup(
    H::AbstractMatrix{<:Integer},
    s::AbstractVector{<:Integer},
    soft_llr::AbstractVector{<:Real},
)
    R, s_reduced, row_pivots = gf2_rref_with_rhs(H, s)
    rank_H = length(row_pivots)
    num_bits = size(H, 2)

    H_reduced = rank_H == 0 ? zeros(Int64, 0, num_bits) : R[1:rank_H, :]
    s_reduced = rank_H == 0 ? Int64[] : s_reduced[1:rank_H]

    order = sortperm(collect(soft_llr); rev=false)

    if rank_H == 0
        basis = Int64[]
    else
        _, ordered_pivots = gf2_rref(H_reduced[:, order])
        basis = Int64.(order[ordered_pivots])
    end

    basis_set = Set(basis)
    remainder = Int64[j for j in order if !(j in basis_set)]

    if rank_H == 0
        base_basis_solution = Int64[]
        delta_basis = zeros(Int64, 0, length(remainder))
    else
        H_basis = H_reduced[:, basis]
        H_basis_inv = gf2_inverse(H_basis)
        base_basis_solution = mod.(H_basis_inv * s_reduced, 2)
        delta_basis = isempty(remainder) ? zeros(Int64, rank_H, 0) : mod.(H_basis_inv * H_reduced[:, remainder], 2)
    end

    return (
        rank=rank_H,
        basis=basis,
        remainder=remainder,
        base_basis_solution=base_basis_solution,
        delta_basis=delta_basis,
        order=order,
    )
end

function _assemble_osd_candidate(
    num_bits::Int,
    basis::Vector{Int64},
    remainder::Vector{Int64},
    basis_bits::AbstractVector{<:Integer},
    remainder_positions::AbstractVector{<:Integer},
)
    candidate = zeros(Int64, num_bits)
    for (i, bit_idx) in enumerate(basis)
        candidate[bit_idx] = basis_bits[i]
    end
    for pos in remainder_positions
        candidate[remainder[pos]] = 1
    end
    return candidate
end

function _for_each_combination(n::Int, k::Int, f::Function)
    if k < 0 || k > n
        return nothing
    end
    if k == 0
        f(Int64[])
        return nothing
    end

    chosen = Vector{Int64}(undef, k)

    function recurse(start_idx::Int, depth::Int)
        if depth > k
            f(copy(chosen))
            return
        end

        stop_idx = n - (k - depth)
        for i in start_idx:stop_idx
            chosen[depth] = i
            recurse(i + 1, depth + 1)
        end
    end

    recurse(1, 1)
    return nothing
end

"""
    osd_decode(H::AbstractMatrix{<:Integer}, s::AbstractVector{<:Integer},
               soft_llr::AbstractVector{<:Real}; osd_order::Int=0, λ::Int=60)

Return the OSD post-processing solution for the binary syndrome equation
`H * e = s (mod 2)` using the BP soft-decision vector `soft_llr`.

This function implements:
- OSD-0 when `osd_order = 0`,
- combination-sweep OSD up to generic order `osd_order` when `osd_order > 0`.

The search rule is:
- order 0: the standard OSD-0 candidate;
- order 1: all singleton flips on the full OSD remainder set;
- orders `2, ..., osd_order`: all combinations on the first `λ` remainder bits.

The return value is a named tuple with fields
- `error`: the final OSD candidate,
- `basis`: the OSD basis set,
- `remainder`: the ordered remainder set,
- `weight`: the Hamming weight of the returned candidate.
"""
function osd_decode(
    H::AbstractMatrix{<:Integer},
    s::AbstractVector{<:Integer},
    soft_llr::AbstractVector{<:Real};
    osd_order::Int=0,
    λ::Int=60,
)
    H = mod.(Int64.(H), 2)
    s = mod.(Int64.(s), 2)

    num_checks, num_bits = size(H)
    if length(s) != num_checks
        error("The syndrome length has to equal the number of rows of H.")
    end
    if length(soft_llr) != num_bits
        error("The soft-decision vector length has to equal the number of columns of H.")
    end
    if osd_order < 0
        error("osd_order has to be nonnegative.")
    end
    if λ < 0
        error("λ has to be nonnegative.")
    end

    setup = _bp_osd_setup(H, s, soft_llr)
    basis = setup.basis
    remainder = setup.remainder
    rank_H = setup.rank
    kprime = length(remainder)

    best_basis_bits = copy(setup.base_basis_solution)
    best_remainder_positions = Int64[]
    best_weight = count(x -> x != 0, setup.base_basis_solution)

    function consider_positions(positions::Vector{Int64})
        candidate_basis_bits = copy(setup.base_basis_solution)
        for pos in positions
            if rank_H > 0
                candidate_basis_bits = mod.(candidate_basis_bits .+ setup.delta_basis[:, pos], 2)
            end
        end

        candidate_weight = count(x -> x != 0, candidate_basis_bits) + length(positions)
        if candidate_weight < best_weight
            best_basis_bits = candidate_basis_bits
            best_remainder_positions = copy(positions)
            best_weight = candidate_weight
        end
        return nothing
    end

    if osd_order >= 1
        for pos in 1:kprime
            consider_positions(Int64[pos])
        end
    end

    if osd_order >= 2
        λeff = min(λ, kprime)
        for order in 2:min(osd_order, λeff)
            _for_each_combination(λeff, order) do positions
                consider_positions(positions)
            end
        end
    end

    candidate = _assemble_osd_candidate(num_bits, basis, remainder, best_basis_bits, best_remainder_positions)

    if mod.(H * candidate, 2) != s
        error("Internal error: the OSD candidate does not satisfy the syndrome equation.")
    end

    return (error=candidate, basis=basis, remainder=remainder, weight=best_weight)
end

"""
    bp_osd_cs_decode(H::AbstractMatrix{<:Integer}, s::AbstractVector{<:Integer}, p::Real;
                     max_iter::Int=size(H, 2),
                     λ::Int=60,
                     osd_order::Int=0,
                     bp_method::Symbol=:ms)

Decode the binary syndrome equation `H * e = s (mod 2)` using the BP+OSD
framework of Roffe, White, Burton, and Campbell.

The decoder first runs a BP stage, then:
- if BP converges to a syndrome-matching solution, that solution is returned;
- otherwise, OSD post-processing is applied.

The OSD stage supports:
- `osd_order = 0`: OSD-0
- `osd_order > 0`: combination-sweep OSD up to order `osd_order`

The BP stage is customizable through `bp_method`:
- `:ms`: scaled min-sum BP (default; this matches the BP variant used by Roffe et al.)

Keyword arguments:
- `max_iter`: maximum number of BP iterations
- `λ`: size cutoff used by the combination-sweep search for orders `≥ 2`
- `osd_order`: OSD order
- `bp_method`: BP variant to use

The return value is a named tuple with fields
- `error`: the final correction
- `converged`: whether BP alone converged
- `bp_error`: the BP hard decision
- `llr`: the final BP soft-decision / reliability vector used by OSD
- `basis`: the OSD basis set
- `remainder`: the ordered OSD remainder set
- `weight`: the Hamming weight of the returned correction
- `bp_method`: the BP method used
"""
function bp_osd_cs_decode(
    H::AbstractMatrix{<:Integer},
    s::AbstractVector{<:Integer},
    p::Real;
    max_iter::Int=size(H, 2),
    λ::Int=60,
    osd_order::Int=0,
    bp_method::Symbol=:ms,
)
    H = mod.(Int64.(H), 2)
    s = mod.(Int64.(s), 2)

    num_checks, num_bits = size(H)
    if length(s) != num_checks
        error("The syndrome length has to equal the number of rows of H.")
    end
    if osd_order < 0
        error("osd_order has to be nonnegative.")
    end
    if λ < 0
        error("λ has to be nonnegative.")
    end

    bp =
        if bp_method == :ms
            bp_min_sum_decode(H, s, p; max_iter=max_iter)
        else
            error("Unsupported bp_method. Check docstring.")
        end

    if bp.converged
        return (
            error=bp.error,
            converged=true,
            bp_error=bp.error,
            llr=bp.llr,
            basis=Int64[],
            remainder=Int64[],
            weight=count(x -> x != 0, bp.error),
            bp_method=bp_method,
        )
    end

    osd = osd_decode(H, s, bp.llr; osd_order=osd_order, λ=λ)

    return (
        error=osd.error,
        converged=false,
        bp_error=bp.error,
        llr=bp.llr,
        basis=osd.basis,
        remainder=osd.remainder,
        weight=osd.weight,
        bp_method=bp_method,
    )
end

"""
    css_bp_osd_cs_decode(HX::AbstractMatrix{<:Integer}, HZ::AbstractMatrix{<:Integer},
                         sx::AbstractVector{<:Integer}, sz::AbstractVector{<:Integer}, p::Real;
                         max_iter::Int=max(size(HX, 2), size(HZ, 2)),
                         λ::Int=60,
                         osd_order::Int=0,
                         bp_method::Symbol=:ms)

Decode a CSS code under uncorrelated code-capacity `X/Z` noise using two
independent calls to `bp_osd_cs_decode`.

The convention is:
- `sx = HZ * x (mod 2)` is the syndrome induced by `X` errors
- `sz = HX * z (mod 2)` is the syndrome induced by `Z` errors

The BP stage is customizable through `bp_method`, see docstring of `bp_osd_cs_decode`

The OSD stage supports:
- `osd_order = 0`: OSD-0
- `osd_order > 0`: combination-sweep OSD up to order `osd_order`

Keyword arguments:
- `max_iter`: maximum number of BP iterations used in each sector
- `λ`: size cutoff used by the combination-sweep search for orders `≥ 2`
- `osd_order`: OSD order used in each sector
- `bp_method`: BP variant used in each sector

The return value is a named tuple with fields
- `x`: the estimated `X`-component correction
- `z`: the estimated `Z`-component correction
- `x_result`: the full decoder output for the `X` sector
- `z_result`: the full decoder output for the `Z` sector
"""
function css_bp_osd_cs_decode(
    HX::AbstractMatrix{<:Integer},
    HZ::AbstractMatrix{<:Integer},
    sx::AbstractVector{<:Integer},
    sz::AbstractVector{<:Integer},
    p::Real;
    max_iter::Int=max(size(HX, 2), size(HZ, 2)),
    λ::Int=60,
    osd_order::Int=0,
    bp_method::Symbol=:ms,
)
    x_result = bp_osd_cs_decode(
        HZ, sx, p;
        max_iter=max_iter,
        λ=λ,
        osd_order=osd_order,
        bp_method=bp_method,
    )

    z_result = bp_osd_cs_decode(
        HX, sz, p;
        max_iter=max_iter,
        λ=λ,
        osd_order=osd_order,
        bp_method=bp_method,
    )

    return (
        x=x_result.error,
        z=z_result.error,
        x_result=x_result,
        z_result=z_result,
    )
end