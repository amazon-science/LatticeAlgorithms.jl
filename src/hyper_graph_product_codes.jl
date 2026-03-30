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

using Random

"""
    hgp_full_rank_repetition_check_matrix(n::Int)

Return the full-rank repetition-code parity-check matrix of size `(n - 1) × n`.

This is the open-chain repetition code with checks
`x₁ + x₂, x₂ + x₃, ..., x_{n-1} + x_n` over `F₂`.
"""
function hgp_full_rank_repetition_check_matrix(n::Int)
    if n < 2
        error("The repetition-code length n has to satisfy n >= 2.")
    end

    H = zeros(Int64, n - 1, n)
    for i in 1:(n - 1)
        H[i, i] = 1
        H[i, i + 1] = 1
    end
    return H
end

"""
    hgp_ring_repetition_check_matrix(n::Int)

Return the closed-loop repetition-code parity-check matrix of size `n × n`.

This is the cycle/ring repetition code with checks
`x₁ + x₂, x₂ + x₃, ..., x_n + x₁` over `F₂`.
"""
function hgp_ring_repetition_check_matrix(n::Int)
    if n < 2
        error("The repetition-code length n has to satisfy n >= 2.")
    end

    H = zeros(Int64, n, n)
    for i in 1:n
        H[i, i] = 1
        H[i, mod1(i + 1, n)] = 1
    end
    return H
end

"""
    hgp_check_matrices(H1::AbstractMatrix{<:Integer}, H2::AbstractMatrix{<:Integer})

Return the CSS check matrices `(HX, HZ)` of the hypergraph-product code built
from the two classical parity-check matrices `H1` and `H2`.

The convention is
- `H1` has size `r1 × n1`,
- `H2` has size `r2 × n2`,
- `HX = [H1 ⊗ I_{n2}   I_{r1} ⊗ H2^T]`,
- `HZ = [I_{n1} ⊗ H2   H1^T ⊗ I_{r2}]`.
"""
function hgp_check_matrices(
    H1::AbstractMatrix{<:Integer},
    H2::AbstractMatrix{<:Integer},
)
    H1 = mod.(Int64.(H1), 2)
    H2 = mod.(Int64.(H2), 2)

    r1, n1 = size(H1)
    r2, n2 = size(H2)

    HX = hcat(
        kron(H1, Matrix{Int64}(I, n2, n2)),
        kron(Matrix{Int64}(I, r1, r1), transpose(H2)),
    )
    HZ = hcat(
        kron(Matrix{Int64}(I, n1, n1), H2),
        kron(transpose(H1), Matrix{Int64}(I, r2, r2)),
    )

    return mod.(Matrix{Int64}(HX), 2), mod.(Matrix{Int64}(HZ), 2)
end

"""
    hgp_num_data_qubits(H1::AbstractMatrix{<:Integer}, H2::AbstractMatrix{<:Integer})

Return the number of physical qubits of the hypergraph-product code.
"""
hgp_num_data_qubits(H1::AbstractMatrix{<:Integer}, H2::AbstractMatrix{<:Integer}) =
    size(H1, 2) * size(H2, 2) + size(H1, 1) * size(H2, 1)

"""
    hgp_num_logicals(H1::AbstractMatrix{<:Integer}, H2::AbstractMatrix{<:Integer})

Return the number of logical qubits of the hypergraph-product code.

This uses the standard dimension formula
`k = k1*k2 + k1T*k2T`, where `k1T` and `k2T` are the dimensions of the
transpose codes.
"""
function hgp_num_logicals(
    H1::AbstractMatrix{<:Integer},
    H2::AbstractMatrix{<:Integer},
)
    k1 = size(gf2_nullspace(mod.(Int64.(H1), 2)), 1)
    k2 = size(gf2_nullspace(mod.(Int64.(H2), 2)), 1)
    k1T = size(gf2_nullspace(mod.(Int64.(transpose(H1)), 2)), 1)
    k2T = size(gf2_nullspace(mod.(Int64.(transpose(H2)), 2)), 1)
    return k1 * k2 + k1T * k2T
end

"""
    hgp_parameters(H1::AbstractMatrix{<:Integer}, H2::AbstractMatrix{<:Integer}; compute_distance=false)

Return a named tuple with the hypergraph-product-code parameters.

By default this function returns `(n, k)`. If `compute_distance=true`, it
returns `(n, k, d)`, where `d` is computed by an exact exhaustive search via
`css_distance(HX, HZ)`.
"""
function hgp_parameters(
    H1::AbstractMatrix{<:Integer},
    H2::AbstractMatrix{<:Integer};
    compute_distance::Bool=false,
)
    n = hgp_num_data_qubits(H1, H2)
    k = hgp_num_logicals(H1, H2)

    if !compute_distance
        return (n=n, k=k)
    end

    HX, HZ = hgp_check_matrices(H1, H2)
    dX = css_distance(HX, HZ; sector=:X)
    dZ = css_distance(HX, HZ; sector=:Z)
    return (n=n, k=k, d=min(dX, dZ))
end

"""
    hgp_X_stabilizers(H1::AbstractMatrix{<:Integer}, H2::AbstractMatrix{<:Integer})

Return the X stabilizers of the hypergraph-product code.
"""
function hgp_X_stabilizers(
    H1::AbstractMatrix{<:Integer},
    H2::AbstractMatrix{<:Integer},
)
    HX, _ = hgp_check_matrices(H1, H2)
    return css_stabilizers_from_check_matrix(HX)
end

"""
    hgp_Z_stabilizers(H1::AbstractMatrix{<:Integer}, H2::AbstractMatrix{<:Integer})

Return the Z stabilizers of the hypergraph-product code.
"""
function hgp_Z_stabilizers(
    H1::AbstractMatrix{<:Integer},
    H2::AbstractMatrix{<:Integer},
)
    _, HZ = hgp_check_matrices(H1, H2)
    return css_stabilizers_from_check_matrix(HZ)
end

"""
    hgp_stabilizers(H1::AbstractMatrix{<:Integer}, H2::AbstractMatrix{<:Integer})

Return the stabilizers of the hypergraph-product code in the package convention:
- X on qubit `j` is encoded as `2j - 1`
- Z on qubit `j` is encoded as `2j`
"""
function hgp_stabilizers(
    H1::AbstractMatrix{<:Integer},
    H2::AbstractMatrix{<:Integer},
)
    X_dict = hgp_X_stabilizers(H1, H2)
    Z_dict = hgp_Z_stabilizers(H1, H2)

    dict = Dict{Int64, Vector{Int64}}()

    for (key, val) in sort(X_dict)
        dict[key] = [2 * item - 1 for item in val]
    end

    for (key, val) in sort(Z_dict)
        dict[length(X_dict) + key] = [2 * item for item in val]
    end

    return dict
end

"""
    hgp_X_logicals(H1::AbstractMatrix{<:Integer}, H2::AbstractMatrix{<:Integer})

Return a dictionary for the logical X operators of the hypergraph-product code.
"""
function hgp_X_logicals(
    H1::AbstractMatrix{<:Integer},
    H2::AbstractMatrix{<:Integer},
)
    HX, HZ = hgp_check_matrices(H1, H2)
    X_basis, _ = css_logicals(HX, HZ)

    dict = Dict{Int64, Vector{Int64}}()
    for i in 1:size(X_basis, 1)
        dict[i] = Int64.(findall(x -> x != 0, vec(X_basis[i, :])))
    end

    return dict
end

"""
    hgp_Z_logicals(H1::AbstractMatrix{<:Integer}, H2::AbstractMatrix{<:Integer})

Return a dictionary for the logical Z operators of the hypergraph-product code.
"""
function hgp_Z_logicals(
    H1::AbstractMatrix{<:Integer},
    H2::AbstractMatrix{<:Integer},
)
    HX, HZ = hgp_check_matrices(H1, H2)
    _, Z_basis = css_logicals(HX, HZ)

    dict = Dict{Int64, Vector{Int64}}()
    for i in 1:size(Z_basis, 1)
        dict[i] = Int64.(findall(x -> x != 0, vec(Z_basis[i, :])))
    end

    return dict
end

"""
    hgp_Mq(H1::AbstractMatrix{<:Integer}, H2::AbstractMatrix{<:Integer})

Return the hypergraph-product-code generator in the `q` subspace.
"""
function hgp_Mq(
    H1::AbstractMatrix{<:Integer},
    H2::AbstractMatrix{<:Integer},
)
    HX, _ = hgp_check_matrices(H1, H2)
    return css_generator_from_check_matrix(HX)
end

"""
    hgp_Mp(H1::AbstractMatrix{<:Integer}, H2::AbstractMatrix{<:Integer})

Return the hypergraph-product-code generator in the `p` subspace.
"""
function hgp_Mp(
    H1::AbstractMatrix{<:Integer},
    H2::AbstractMatrix{<:Integer},
)
    _, HZ = hgp_check_matrices(H1, H2)
    return css_generator_from_check_matrix(HZ)
end

"""
    hgp_M(H1::AbstractMatrix{<:Integer}, H2::AbstractMatrix{<:Integer})

Return the full hypergraph-product-code GKP generator matrix.
"""
function hgp_M(
    H1::AbstractMatrix{<:Integer},
    H2::AbstractMatrix{<:Integer},
)
    Mq = hgp_Mq(H1, H2)
    Mp = hgp_Mp(H1, H2)
    M = BlockDiagonal([Mq, Mp])
    T = basis_transformation(size(Mq, 1))
    return transpose(T) * M * T
end


"""
    panteleev_kalachev_c2_parent_check_matrix()

Return the `31 × 31` binary circulant parent check matrix used for the fixed
hypergraph-product code `C2 = [[1922, 50, 16]]` in Panteleev–Kalachev.

The parent cyclic code has length `ℓ = 31` and parity polynomial
`h(x) = 1 + x^2 + x^5`.
"""
panteleev_kalachev_c2_parent_check_matrix() = circulant_binary_matrix(31, [0, 2, 5])
