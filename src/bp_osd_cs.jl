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
    _bp_osd_setup(H::AbstractMatrix{<:Integer},
                  s::AbstractVector{<:Integer},
                  soft_llr::AbstractVector{<:Real})

Construct the linear-algebra data used by OSD from the BP soft output.

The columns are ordered from least reliable to most reliable using increasing
`abs.(soft_llr)`. The returned named tuple contains

- `rank`
- `basis`
- `remainder`
- `base_basis_solution`
- `delta_basis`
- `order`
"""
function _bp_osd_setup(
    H::AbstractMatrix{<:Integer},
    s::AbstractVector{<:Integer},
    soft_llr::AbstractVector{<:Real},
)
    H = mod.(Int64.(H), 2)
    s = mod.(Int64.(s), 2)
    soft_llr = Float64.(collect(soft_llr))

    num_checks, num_bits = size(H)
    if length(s) != num_checks
        error("The syndrome length has to equal the number of rows of H.")
    end
    if length(soft_llr) != num_bits
        error("The soft-decision vector length has to equal the number of columns of H.")
    end

    R, s_reduced, row_pivots = gf2_rref_with_rhs(H, s)
    rank_H = length(row_pivots)

    H_reduced = rank_H == 0 ? zeros(Int64, 0, num_bits) : R[1:rank_H, :]
    s_reduced = rank_H == 0 ? Int64[] : s_reduced[1:rank_H]

    # Reliability ordering: small |LLR| = less reliable.
    order = sortperm(abs.(soft_llr); rev=false)

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
        delta_basis = isempty(remainder) ? zeros(Int64, rank_H, 0) :
            mod.(H_basis_inv * H_reduced[:, remainder], 2)
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

"""
    _for_each_combination(n::Int, k::Int, f::Function)

Call `f` on every `k`-subset of `1:n`, represented as a vector of positions.
"""
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
    osd_decode(H::AbstractMatrix{<:Integer},
               s::AbstractVector{<:Integer},
               soft_llr::AbstractVector{<:Real};
               osd_order::Int=0,
               λ::Int=60)

Run ordered-statistics decoding (OSD) using the supplied BP soft output.

This function implements
- OSD-0 when `osd_order = 0`,
- combination-sweep OSD up to generic order `osd_order` when `osd_order > 0`.

The search rule is
- order 0: the standard OSD-0 candidate;
- order 1: all singleton flips on the full OSD remainder set;
- orders `2, ..., osd_order`: all combinations on the first `λ` remainder bits.

Return value
------------
A named tuple with fields
- `error`
- `basis`
- `remainder`
- `remainder_positions`
- `hamming_weight`
- `reliability_order`
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
    best_hamming_weight = count(!iszero, setup.base_basis_solution)

    function consider_positions(positions::Vector{Int64})
        candidate_basis_bits = copy(setup.base_basis_solution)

        for pos in positions
            if rank_H > 0
                candidate_basis_bits = mod.(candidate_basis_bits .+ setup.delta_basis[:, pos], 2)
            end
        end

        candidate_hamming_weight =
            count(!iszero, candidate_basis_bits) + length(positions)

        if candidate_hamming_weight < best_hamming_weight
            best_basis_bits = candidate_basis_bits
            best_remainder_positions = copy(positions)
            best_hamming_weight = candidate_hamming_weight
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

    # Assemble the full-length OSD candidate from its basis coordinates and the chosen
    # positions in the remainder set
    candidate = zeros(Int64, num_bits)
    for (i, bit_idx) in enumerate(basis)
        candidate[bit_idx] = best_basis_bits[i]
    end
    for pos in best_remainder_positions
        candidate[remainder[pos]] = 1
    end
    return candidate
end    

    if mod.(H * candidate, 2) != s
        error("Internal error: the OSD candidate does not satisfy the syndrome equation.")
    end

    return (
        error=candidate,
        basis=basis,
        remainder=remainder,
        remainder_positions=best_remainder_positions,
        hamming_weight=best_hamming_weight,
        reliability_order=setup.order,
    )
end

"""
    bp_osd_cs_decode(H::AbstractMatrix{<:Integer},
                     s::AbstractVector{<:Integer},
                     p;
                     max_iter::Int=size(H, 2),
                     λ::Int=60,
                     osd_order::Int=0,
                     check_to_bit_update_rule::Symbol=:min_sum,
                     bit_to_check_update_rule::Symbol=:memoryless,
                     γ=nothing,
                     initial_marginals=nothing,
                     min_sum_scaling::Symbol=:roffe)

Decode the binary syndrome equation `H * e = s (mod 2)` using BP followed by
OSD post-processing if BP alone does not converge.

This function is a thin BP+OSD wrapper around `bp_decode` and `osd_decode`.
Its BP keywords intentionally mirror those of `bp_decode`. The default BP
settings match the min-sum / memoryless / Roffe-style scaling path traditionally
used in BP+OSD.

The decoder proceeds as follows:
1. run `bp_decode`;
2. if BP returns a syndrome-consistent hard decision, return it immediately;
3. otherwise, run `osd_decode` on the final BP `llr`.

Return value
------------
A named tuple with fields
- `error`: final returned correction
- `converged`: whether BP alone converged
- `bp_result`: full result returned by `bp_decode`
- `osd_result`: full result returned by `osd_decode`, or `nothing`
- `llr`: BP soft output used by OSD
- `basis`
- `remainder`
- `hamming_weight`
"""
function bp_osd_cs_decode(
    H::AbstractMatrix{<:Integer},
    s::AbstractVector{<:Integer},
    p;
    max_iter::Int=size(H, 2),
    λ::Int=60,
    osd_order::Int=0,
    check_to_bit_update_rule::Symbol=:min_sum,
    bit_to_check_update_rule::Symbol=:memoryless,
    γ=nothing,
    initial_marginals=nothing,
    min_sum_scaling::Symbol=:roffe,
)
    H = mod.(Int64.(H), 2)
    s = mod.(Int64.(s), 2)

    num_checks, _ = size(H)
    if length(s) != num_checks
        error("The syndrome length has to equal the number of rows of H.")
    end
    if osd_order < 0
        error("osd_order has to be nonnegative.")
    end
    if λ < 0
        error("λ has to be nonnegative.")
    end

    bp = bp_decode(
        H,
        s,
        p;
        max_iter=max_iter,
        check_to_bit_update_rule=check_to_bit_update_rule,
        bit_to_check_update_rule=bit_to_check_update_rule,
        γ=γ,
        initial_marginals=initial_marginals,
        min_sum_scaling=min_sum_scaling,
    )

    if bp.converged
        return (
            error=bp.error,
            converged=true,
            bp_result=bp,
            osd_result=nothing,
            llr=bp.llr,
            basis=Int64[],
            remainder=Int64[],
            hamming_weight=count(!iszero, bp.error),
        )
    end

    osd = osd_decode(H, s, bp.llr; osd_order=osd_order, λ=λ)

    return (
        error=osd.error,
        converged=false,
        bp_result=bp,
        osd_result=osd,
        llr=bp.llr,
        basis=osd.basis,
        remainder=osd.remainder,
        hamming_weight=osd.hamming_weight,
    )
end

"""
    css_bp_osd_cs_decode(HX::AbstractMatrix{<:Integer},
                         HZ::AbstractMatrix{<:Integer},
                         sx::AbstractVector{<:Integer},
                         sz::AbstractVector{<:Integer},
                         px,
                         pz=px;
                         max_iter::Int=max(size(HX, 2), size(HZ, 2)),
                         λ::Int=60,
                         osd_order::Int=0,
                         check_to_bit_update_rule::Symbol=:min_sum,
                         bit_to_check_update_rule::Symbol=:memoryless,
                         γx=nothing,
                         γz=nothing,
                         initial_marginals_x=nothing,
                         initial_marginals_z=nothing,
                         min_sum_scaling::Symbol=:roffe)

Decode a CSS code under independent code-capacity `X/Z` noise using two
independent calls to `bp_osd_cs_decode`.

The convention is:
- `sx = HZ * x (mod 2)` is the syndrome induced by `X` errors,
- `sz = HX * z (mod 2)` is the syndrome induced by `Z` errors.

`px` and `pz` can each be either a scalar bit-flip probability or a vector of
per-bit probabilities, matching the interface of `bp_decode`.

Return value
------------
A named tuple with fields
- `x`
- `z`
- `x_result`
- `z_result`
"""
function css_bp_osd_cs_decode(
    HX::AbstractMatrix{<:Integer},
    HZ::AbstractMatrix{<:Integer},
    sx::AbstractVector{<:Integer},
    sz::AbstractVector{<:Integer},
    px,
    pz=px;
    max_iter::Int=max(size(HX, 2), size(HZ, 2)),
    λ::Int=60,
    osd_order::Int=0,
    check_to_bit_update_rule::Symbol=:min_sum,
    bit_to_check_update_rule::Symbol=:memoryless,
    γx=nothing,
    γz=nothing,
    initial_marginals_x=nothing,
    initial_marginals_z=nothing,
    min_sum_scaling::Symbol=:roffe,
)
    x_result = bp_osd_cs_decode(
        HZ,
        sx,
        px;
        max_iter=max_iter,
        λ=λ,
        osd_order=osd_order,
        check_to_bit_update_rule=check_to_bit_update_rule,
        bit_to_check_update_rule=bit_to_check_update_rule,
        γ=γx,
        initial_marginals=initial_marginals_x,
        min_sum_scaling=min_sum_scaling,
    )

    z_result = bp_osd_cs_decode(
        HX,
        sz,
        pz;
        max_iter=max_iter,
        λ=λ,
        osd_order=osd_order,
        check_to_bit_update_rule=check_to_bit_update_rule,
        bit_to_check_update_rule=bit_to_check_update_rule,
        γ=γz,
        initial_marginals=initial_marginals_z,
        min_sum_scaling=min_sum_scaling,
    )

    return (
        x=x_result.error,
        z=z_result.error,
        x_result=x_result,
        z_result=z_result,
    )
end