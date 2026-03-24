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
    gf2_rref(M::AbstractMatrix{<:Integer})

Return the reduced row-echelon form of `M` over `F_2`, together with the pivot
columns.

Output:
- `R`: the reduced row-echelon form of `M` over `F_2`
- `pivots`: the vector of pivot columns

Example:
```julia
R, pivots = gf2_rref([1 1 0; 0 1 1])
# R = [1 0 1; 0 1 1]
# pivots = [1, 2]
```
"""
function gf2_rref(M::AbstractMatrix{<:Integer})
    R = mod.(Matrix{Int64}(M), 2)
    num_rows, num_cols = size(R)

    pivots = Int64[]
    row_idx = 1

    for col_idx in 1:num_cols
        pivot_row = nothing
        for i in row_idx:num_rows
            if R[i, col_idx] == 1
                pivot_row = i
                break
            end
        end

        if isnothing(pivot_row)
            continue
        end

        if pivot_row != row_idx
            R[row_idx, :], R[pivot_row, :] = copy(R[pivot_row, :]), copy(R[row_idx, :])
        end

        for i in 1:num_rows
            if i != row_idx && R[i, col_idx] == 1
                R[i, :] = mod.(R[i, :] .+ R[row_idx, :], 2)
            end
        end

        push!(pivots, col_idx)
        row_idx += 1

        if row_idx > num_rows
            break
        end
    end

    return R, pivots
end

"""
    gf2_nullspace(M::AbstractMatrix{<:Integer})

Return a row-basis for the nullspace of `M` over `F_2`.

If the nullspace is trivial, this function returns a `0 × n` matrix, where
`n = size(M, 2)`.

Example:
```julia
N = gf2_nullspace([1 1 0; 0 1 1])
# N = [1 1 1]
```
"""
function gf2_nullspace(M::AbstractMatrix{<:Integer})
    R, pivots = gf2_rref(M)
    num_cols = size(M, 2)
    free_cols = [j for j in 1:num_cols if !(j in pivots)]

    if isempty(free_cols)
        return zeros(Int64, 0, num_cols)
    end

    basis = zeros(Int64, length(free_cols), num_cols)

    for (i, free_col) in enumerate(free_cols)
        basis[i, free_col] = 1
        for (row_idx, pivot_col) in enumerate(pivots)
            basis[i, pivot_col] = R[row_idx, free_col]
        end
    end

    return basis
end

"""
    gf2_inverse(M::AbstractMatrix{<:Integer})

Return the inverse of the square matrix `M` over `F_2`.

An error is thrown if `M` is not square or is singular over `F_2`.

Example:
```julia
Minv = gf2_inverse([1 1; 1 0])
@test mod.([1 1; 1 0] * Minv, 2) == Matrix{Int64}(I, 2, 2)
```
"""
function gf2_inverse(M::AbstractMatrix{<:Integer})
    num_rows, num_cols = size(M)

    if num_rows != num_cols
        error("The input matrix has to be square.")
    end

    A = hcat(mod.(Matrix{Int64}(M), 2), Matrix{Int64}(I, num_rows, num_rows))

    row_idx = 1
    for col_idx in 1:num_cols
        pivot_row = nothing
        for i in row_idx:num_rows
            if A[i, col_idx] == 1
                pivot_row = i
                break
            end
        end

        if isnothing(pivot_row)
            error("The matrix is singular over F_2.")
        end

        if pivot_row != row_idx
            A[row_idx, :], A[pivot_row, :] = copy(A[pivot_row, :]), copy(A[row_idx, :])
        end

        for i in 1:num_rows
            if i != row_idx && A[i, col_idx] == 1
                A[i, :] = mod.(A[i, :] .+ A[row_idx, :], 2)
            end
        end

        row_idx += 1
    end

    return A[:, num_cols+1:end]
end

"""
    css_stabilizers_from_check_matrix(H::AbstractMatrix{<:Integer})

Return the stabilizer-support dictionary associated with the binary check matrix
`H`.

The `i`-th row of `H` is converted into the support of the `i`-th stabilizer,
stored as a vector of 1-based qubit indices.

Example:
```julia
css_stabilizers_from_check_matrix([1 0 1; 0 1 1])
# Dict(1 => [1, 3], 2 => [2, 3])
```
"""
function css_stabilizers_from_check_matrix(H::AbstractMatrix{<:Integer})
    dict = Dict{Int64, Vector{Int64}}()
    for row_idx in 1:size(H, 1)
        dict[row_idx] = Int64.(findall(x -> x != 0, vec(H[row_idx, :])))
    end
    return dict
end

"""
    css_generator_from_check_matrix(H::AbstractMatrix{<:Integer})

Return the GKP generator in one quadrature associated with the CSS check matrix
`H`.

The rows of `H` are reduced to an independent row-basis over `F_2`, and this
basis is embedded into a matrix `M` such that `det(M) = √2` after dividing by
`√2`.

This is the generic utility used by `Mq` and `Mp` constructors for CSS LDPC
families.
"""
function css_generator_from_check_matrix(H::AbstractMatrix{<:Integer})
    R, pivots = gf2_rref(H)
    num_qubits = size(H, 2)

    M = 2 * Matrix{Int64}(I, num_qubits, num_qubits)

    for (row_idx, pivot_col) in enumerate(pivots)
        M[pivot_col, :] = R[row_idx, :]
    end

    return M / √2
end

"""
    css_logicals(HX::AbstractMatrix{<:Integer}, HZ::AbstractMatrix{<:Integer})

Return paired row-bases `(X_basis, Z_basis)` for the logical `X` and `Z`
operators of a CSS code specified by the binary check matrices `HX` and `HZ`.

The returned matrices satisfy:
- each row of `X_basis` commutes with all rows of `HZ`
- each row of `Z_basis` commutes with all rows of `HX`
- the pairing `X_basis * transpose(Z_basis)` is the identity over `F_2`

Each row of the output matrices is a binary indicator vector on the physical
qubits.
"""
function css_logicals(HX::AbstractMatrix{<:Integer}, HZ::AbstractMatrix{<:Integer})
    num_qubits = size(HX, 2)

    RX, pivots_X = gf2_rref(HX)
    RZ, pivots_Z = gf2_rref(HZ)

    rowspace_X = isempty(pivots_X) ? zeros(Int64, 0, num_qubits) : RX[1:length(pivots_X), :]
    rowspace_Z = isempty(pivots_Z) ? zeros(Int64, 0, num_qubits) : RZ[1:length(pivots_Z), :]

    nullspace_HX = gf2_nullspace(HX)
    nullspace_HZ = gf2_nullspace(HZ)

    basis_Z = copy(rowspace_Z)
    Z_basis = zeros(Int64, 0, num_qubits)
    rank_Z = size(basis_Z, 1)

    for i in 1:size(nullspace_HX, 1)
        candidate = reshape(nullspace_HX[i, :], 1, :)
        _, pivots = gf2_rref(vcat(basis_Z, candidate))
        if length(pivots) > rank_Z
            basis_Z = vcat(basis_Z, candidate)
            Z_basis = vcat(Z_basis, candidate)
            rank_Z = length(pivots)
        end
    end

    basis_X = copy(rowspace_X)
    X_basis = zeros(Int64, 0, num_qubits)
    rank_X = size(basis_X, 1)

    for i in 1:size(nullspace_HZ, 1)
        candidate = reshape(nullspace_HZ[i, :], 1, :)
        _, pivots = gf2_rref(vcat(basis_X, candidate))
        if length(pivots) > rank_X
            basis_X = vcat(basis_X, candidate)
            X_basis = vcat(X_basis, candidate)
            rank_X = length(pivots)
        end
    end

    if size(X_basis, 1) != size(Z_basis, 1)
        error("Inconsistent number of X and Z logical operators.")
    end

    if size(X_basis, 1) == 0
        return X_basis, Z_basis
    end

    pairing = mod.(X_basis * transpose(Z_basis), 2)
    X_basis = mod.(gf2_inverse(pairing) * X_basis, 2)

    return X_basis, Z_basis
end
