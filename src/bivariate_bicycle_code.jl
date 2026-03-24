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
    bb_cyclic_shift_matrix(n::Int)

Return the cyclic shift matrix `S_n` of size `n × n` over `F_2`.
"""
function bb_cyclic_shift_matrix(n::Int)
    if n < 1
        error("The matrix size n has to be a positive integer.")
    end

    S = zeros(Int64, n, n)
    for i in 1:n
        S[i, mod1(i + 1, n)] = 1
    end

    return S
end

"""
    bb_x_matrix(l::Int, m::Int)

Return the matrix `x = S_l ⊗ I_m` used in the bivariate bicycle code
construction.
"""
function bb_x_matrix(l::Int, m::Int)
    if l < 1 || m < 1
        error("Both l and m have to be positive integers.")
    end
    return kron(bb_cyclic_shift_matrix(l), Matrix{Int64}(I, m, m))
end

"""
    bb_y_matrix(l::Int, m::Int)

Return the matrix `y = I_l ⊗ S_m` used in the bivariate bicycle code
construction.
"""
function bb_y_matrix(l::Int, m::Int)
    if l < 1 || m < 1
        error("Both l and m have to be positive integers.")
    end
    return kron(Matrix{Int64}(I, l, l), bb_cyclic_shift_matrix(m))
end

function _bb_matrix(l::Int, m::Int, terms::Vector{Tuple{Symbol, Int}})
    if length(terms) != 3
        error("A BB polynomial must contain exactly three monomials.")
    end

    normalized_terms = Tuple{Symbol, Int64}[]
    for (kind, power) in terms
        if kind == :I
            push!(normalized_terms, (:I, 0))
        elseif kind == :x
            p = mod(power, l)
            push!(normalized_terms, p == 0 ? (:I, 0) : (:x, p))
        elseif kind == :y
            p = mod(power, m)
            push!(normalized_terms, p == 0 ? (:I, 0) : (:y, p))
        else
            error("Each BB monomial has to be of the form (:I, 0), (:x, p), or (:y, p).")
        end
    end

    if length(unique(normalized_terms)) != 3
        error("The monomials are not distinct after reducing exponents modulo l or m.")
    end

    x = bb_x_matrix(l, m)
    y = bb_y_matrix(l, m)
    M = zeros(Int64, l * m, l * m)

    for (kind, power) in normalized_terms
        if kind == :I
            M .+= Matrix{Int64}(I, l * m, l * m)
        elseif kind == :x
            M .+= x^power
        else
            M .+= y^power
        end
    end

    return mod.(M, 2)
end

"""
    bb_matrices(l::Int, m::Int, A_terms::Vector{Tuple{Symbol, Int}}, B_terms::Vector{Tuple{Symbol, Int}})

Return the pair `(A, B)` defining the BB code.
"""
function bb_matrices(
    l::Int,
    m::Int,
    A_terms::Vector{Tuple{Symbol, Int}},
    B_terms::Vector{Tuple{Symbol, Int}}
)
    A = _bb_matrix(l, m, A_terms)
    B = _bb_matrix(l, m, B_terms)
    return A, B
end

"""
    bb_check_matrices(l::Int, m::Int, A_terms::Vector{Tuple{Symbol, Int}}, B_terms::Vector{Tuple{Symbol, Int}})

Return the BB CSS check matrices `(HX, HZ)`.
"""
function bb_check_matrices(
    l::Int,
    m::Int,
    A_terms::Vector{Tuple{Symbol, Int}},
    B_terms::Vector{Tuple{Symbol, Int}}
)
    A, B = bb_matrices(l, m, A_terms, B_terms)
    HX = hcat(A, B)
    HZ = hcat(transpose(B), transpose(A))
    return Matrix{Int64}(HX), Matrix{Int64}(HZ)
end

"""
    bb_num_data_qubits(l::Int, m::Int)

Return the number of physical qubits of the BB code.
"""
bb_num_data_qubits(l::Int, m::Int) = 2 * l * m

"""
    bb_num_logicals(l::Int, m::Int, A_terms::Vector{Tuple{Symbol, Int}}, B_terms::Vector{Tuple{Symbol, Int}})

Return the number of logical qubits of the BB code.
"""
function bb_num_logicals(
    l::Int,
    m::Int,
    A_terms::Vector{Tuple{Symbol, Int}},
    B_terms::Vector{Tuple{Symbol, Int}}
)
    A, B = bb_matrices(l, m, A_terms, B_terms)
    return 2 * size(gf2_nullspace(vcat(A, B)), 1)
end

"""
    bb_parameters(l::Int, m::Int, A_terms::Vector{Tuple{Symbol, Int}}, B_terms::Vector{Tuple{Symbol, Int}}; compute_distance=false)

Return a named tuple with the BB-code parameters.

By default this function returns `(n, k)`. If `compute_distance=true`, it
returns `(n, k, d)`, where `d` is computed by an exact exhaustive search via
`css_distance(HX, HZ)`.
"""
function bb_parameters(
    l::Int,
    m::Int,
    A_terms::Vector{Tuple{Symbol, Int}},
    B_terms::Vector{Tuple{Symbol, Int}};
    compute_distance::Bool=false,
)
    n = bb_num_data_qubits(l, m)
    k = bb_num_logicals(l, m, A_terms, B_terms)

    if !compute_distance
        return (n=n, k=k)
    end

    HX, HZ = bb_check_matrices(l, m, A_terms, B_terms)
    d = css_distance(HX, HZ; sector=:Z)
    return (n=n, k=k, d=d)
end

"""
    bb_X_stabilizers(l::Int, m::Int, A_terms::Vector{Tuple{Symbol, Int}}, B_terms::Vector{Tuple{Symbol, Int}})

Return the X stabilizers of the BB code.
"""
function bb_X_stabilizers(
    l::Int,
    m::Int,
    A_terms::Vector{Tuple{Symbol, Int}},
    B_terms::Vector{Tuple{Symbol, Int}}
)
    HX, _ = bb_check_matrices(l, m, A_terms, B_terms)
    return css_stabilizers_from_check_matrix(HX)
end

"""
    bb_Z_stabilizers(l::Int, m::Int, A_terms::Vector{Tuple{Symbol, Int}}, B_terms::Vector{Tuple{Symbol, Int}})

Return the Z stabilizers of the BB code.
"""
function bb_Z_stabilizers(
    l::Int,
    m::Int,
    A_terms::Vector{Tuple{Symbol, Int}},
    B_terms::Vector{Tuple{Symbol, Int}}
)
    _, HZ = bb_check_matrices(l, m, A_terms, B_terms)
    return css_stabilizers_from_check_matrix(HZ)
end

"""
    bb_stabilizers(l::Int, m::Int, A_terms::Vector{Tuple{Symbol, Int}}, B_terms::Vector{Tuple{Symbol, Int}})

Return the stabilizers of the BB code in the package convention:
- X on qubit `j` is encoded as `2j - 1`
- Z on qubit `j` is encoded as `2j`
"""
function bb_stabilizers(
    l::Int,
    m::Int,
    A_terms::Vector{Tuple{Symbol, Int}},
    B_terms::Vector{Tuple{Symbol, Int}}
)
    X_dict = bb_X_stabilizers(l, m, A_terms, B_terms)
    Z_dict = bb_Z_stabilizers(l, m, A_terms, B_terms)

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
    bb_X_logicals(l::Int, m::Int, A_terms::Vector{Tuple{Symbol, Int}}, B_terms::Vector{Tuple{Symbol, Int}})

Return a dictionary for the logical X operators of the BB code.
"""
function bb_X_logicals(
    l::Int,
    m::Int,
    A_terms::Vector{Tuple{Symbol, Int}},
    B_terms::Vector{Tuple{Symbol, Int}}
)
    HX, HZ = bb_check_matrices(l, m, A_terms, B_terms)
    X_basis, _ = css_logicals(HX, HZ)

    dict = Dict{Int64, Vector{Int64}}()
    for i in 1:size(X_basis, 1)
        dict[i] = Int64.(findall(x -> x != 0, vec(X_basis[i, :])))
    end

    return dict
end

"""
    bb_Z_logicals(l::Int, m::Int, A_terms::Vector{Tuple{Symbol, Int}}, B_terms::Vector{Tuple{Symbol, Int}})

Return a dictionary for the logical Z operators of the BB code.
"""
function bb_Z_logicals(
    l::Int,
    m::Int,
    A_terms::Vector{Tuple{Symbol, Int}},
    B_terms::Vector{Tuple{Symbol, Int}}
)
    HX, HZ = bb_check_matrices(l, m, A_terms, B_terms)
    _, Z_basis = css_logicals(HX, HZ)

    dict = Dict{Int64, Vector{Int64}}()
    for i in 1:size(Z_basis, 1)
        dict[i] = Int64.(findall(x -> x != 0, vec(Z_basis[i, :])))
    end

    return dict
end

"""
    bb_Mq(l::Int, m::Int, A_terms::Vector{Tuple{Symbol, Int}}, B_terms::Vector{Tuple{Symbol, Int}})

Return the BB-code generator in the `q` subspace.
"""
function bb_Mq(
    l::Int,
    m::Int,
    A_terms::Vector{Tuple{Symbol, Int}},
    B_terms::Vector{Tuple{Symbol, Int}}
)
    HX, _ = bb_check_matrices(l, m, A_terms, B_terms)
    return css_generator_from_check_matrix(HX)
end

"""
    bb_Mp(l::Int, m::Int, A_terms::Vector{Tuple{Symbol, Int}}, B_terms::Vector{Tuple{Symbol, Int}})

Return the BB-code generator in the `p` subspace.
"""
function bb_Mp(
    l::Int,
    m::Int,
    A_terms::Vector{Tuple{Symbol, Int}},
    B_terms::Vector{Tuple{Symbol, Int}}
)
    _, HZ = bb_check_matrices(l, m, A_terms, B_terms)
    return css_generator_from_check_matrix(HZ)
end

"""
    bb_M(l::Int, m::Int, A_terms::Vector{Tuple{Symbol, Int}}, B_terms::Vector{Tuple{Symbol, Int}})

Return the full BB-code GKP generator matrix.
"""
function bb_M(
    l::Int,
    m::Int,
    A_terms::Vector{Tuple{Symbol, Int}},
    B_terms::Vector{Tuple{Symbol, Int}}
)
    Mq = bb_Mq(l, m, A_terms, B_terms)
    Mp = bb_Mp(l, m, A_terms, B_terms)
    M = BlockDiagonal([Mq, Mp])
    T = basis_transformation(size(Mq, 1))
    return transpose(T) * M * T
end
