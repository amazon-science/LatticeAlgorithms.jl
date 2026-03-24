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


function gf2_rref_with_rhs(H::AbstractMatrix{<:Integer}, s::AbstractVector{<:Integer})
    A = hcat(mod.(Int64.(H), 2), reshape(mod.(Int64.(s), 2), :, 1))
    num_rows, num_cols_plus_rhs = size(A)
    num_cols = num_cols_plus_rhs - 1

    pivots = Int64[]
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
            continue
        end

        if pivot_row != row_idx
            A[row_idx, :], A[pivot_row, :] = copy(A[pivot_row, :]), copy(A[row_idx, :])
        end

        for i in 1:num_rows
            if i != row_idx && A[i, col_idx] == 1
                A[i, :] = mod.(A[i, :] .+ A[row_idx, :], 2)
            end
        end

        push!(pivots, col_idx)
        row_idx += 1

        if row_idx > num_rows
            break
        end
    end

    return A[:, 1:num_cols], vec(A[:, num_cols + 1]), pivots
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

"""
    css_distance(HX::AbstractMatrix{<:Integer}, HZ::AbstractMatrix{<:Integer};
                 max_weight=nothing, sector=:both)

Return the CSS code distance determined by the binary check matrices `HX` and `HZ`.

The keyword `sector` controls which logical sector is searched:
- `sector = :both`: compute both `dX` and `dZ`, and return `min(dX, dZ)`
- `sector = :Z`: compute only the minimum weight of a nontrivial `Z` logical in
  `ker(HX) \\ rowspace(HZ)`
- `sector = :X`: compute only the minimum weight of a nontrivial `X` logical in
  `ker(HZ) \\ rowspace(HX)`

This implementation searches by increasing logical-basis coefficient weight, so
it is exact, but it is only practical when the true distance is reasonably
small. The optional keyword `max_weight` can be used to stop the search early.
In that case:
- if a logical operator of weight at most `max_weight` is found, its exact
  weight is returned;
- otherwise `nothing` is returned.
"""
function css_distance(
    HX::AbstractMatrix{<:Integer},
    HZ::AbstractMatrix{<:Integer};
    max_weight=nothing,
    sector::Symbol=:both,
)
    if size(HX, 2) != size(HZ, 2)
        error("HX and HZ must have the same number of columns.")
    end

    if any(mod.(HX * transpose(HZ), 2) .!= 0)
        error("HX and HZ do not define a CSS code because HX * transpose(HZ) != 0 over F_2.")
    end

    if !(sector in (:both, :Z, :X))
        error("sector must be one of :both, :Z, or :X.")
    end

    num_qubits = size(HX, 2)

    function in_rowspace_rref(v::Vector{Int64}, R::Matrix{Int64}, pivots::Vector{Int64})
        w = copy(v)
        for (row_idx, pivot_col) in enumerate(pivots)
            if w[pivot_col] == 1
                @inbounds for j in 1:length(w)
                    w[j] = xor(w[j], R[row_idx, j])
                end
            end
        end
        return all(x -> x == 0, w)
    end

    function toggle_row!(v::Vector{Int64}, row::AbstractVector{<:Integer})
        @inbounds for j in 1:length(v)
            v[j] = xor(v[j], row[j])
        end
        return nothing
    end

    function min_nontrivial_logical_weight(
        H_kernel::AbstractMatrix{<:Integer},
        H_trivial::AbstractMatrix{<:Integer},
    )
        R_kernel, pivots_kernel = gf2_rref(H_kernel)
        free_cols = Int64[j for j in 1:num_qubits if !(j in pivots_kernel)]

        null_basis = gf2_nullspace(H_kernel)
        num_basis_vectors = size(null_basis, 1)

        if num_basis_vectors == 0
            return nothing
        end

        trivial_coord_rref, trivial_coord_pivots = gf2_rref(H_trivial[:, free_cols])

        best = Ref(isnothing(max_weight) ? (num_qubits + 1) : (max_weight + 1))
        current_codeword = zeros(Int64, num_qubits)
        current_coeffs = zeros(Int64, num_basis_vectors)

        function recurse(start_idx::Int, num_left::Int)
            if num_left == 0
                if !in_rowspace_rref(current_coeffs, trivial_coord_rref, trivial_coord_pivots)
                    weight = count(x -> x != 0, current_codeword)
                    if weight < best[]
                        best[] = weight
                    end
                end
                return
            end

            stop_idx = num_basis_vectors - num_left + 1
            for i in start_idx:stop_idx
                current_coeffs[i] = 1
                toggle_row!(current_codeword, view(null_basis, i, :))
                recurse(i + 1, num_left - 1)
                toggle_row!(current_codeword, view(null_basis, i, :))
                current_coeffs[i] = 0
            end
        end

        coeff_weight = 1
        coeff_weight_limit = isnothing(max_weight) ? num_basis_vectors : min(num_basis_vectors, max_weight)

        while coeff_weight <= coeff_weight_limit && coeff_weight < best[]
            recurse(1, coeff_weight)
            coeff_weight += 1
        end

        return best[] == (isnothing(max_weight) ? num_qubits + 1 : max_weight + 1) ? nothing : best[]
    end

    if sector == :Z
        return min_nontrivial_logical_weight(HX, HZ)
    elseif sector == :X
        return min_nontrivial_logical_weight(HZ, HX)
    else
        dZ = min_nontrivial_logical_weight(HX, HZ)
        dX = min_nontrivial_logical_weight(HZ, HX)

        if isnothing(dX) && isnothing(dZ)
            return nothing
        elseif isnothing(dX)
            return dZ
        elseif isnothing(dZ)
            return dX
        else
            return min(dX, dZ)
        end
    end
end
