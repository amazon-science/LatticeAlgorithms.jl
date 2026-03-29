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

using LatticeAlgorithms
using Test
using LinearAlgebra
using SparseArrays

@testset "css_ldpc_utils" begin
    M = [1 1 0; 0 1 1]
    R, pivots = gf2_rref(M)
    @test R == [1 0 1; 0 1 1]
    @test pivots == [1, 2]

    N = gf2_nullspace(M)
    @test size(N) == (1, 3)
    @test N == [1 1 1]
    @test mod.(M * transpose(N), 2) == zeros(Int64, 2, 1)

    A = [1 1; 1 0]
    Ainv = gf2_inverse(A)
    @test mod.(A * Ainv, 2) == Matrix{Int64}(I, 2, 2)
    @test mod.(Ainv * A, 2) == Matrix{Int64}(I, 2, 2)

    H = [1 0 1; 0 1 1]
    stabs = css_stabilizers_from_check_matrix(H)
    @test stabs == Dict(1 => [1, 3], 2 => [2, 3])

    G = css_generator_from_check_matrix(H)
    @test size(G) == (3, 3)
    @test isapprox(√2 * G, [1 0 1; 0 1 1; 0 0 2])

    HX = [1 1 0; 0 1 1]
    HZ = zeros(Int64, 0, 3)

    X_basis, Z_basis = css_logicals(HX, HZ)

    @test size(X_basis) == (1, 3)
    @test size(Z_basis) == (1, 3)

    @test mod.(HX * transpose(Z_basis), 2) == zeros(Int64, size(HX, 1), size(Z_basis, 1))
    @test mod.(HZ * transpose(X_basis), 2) == zeros(Int64, size(HZ, 1), size(X_basis, 1))
    @test mod.(X_basis * transpose(Z_basis), 2) == Matrix{Int64}(I, 1, 1)
    @test css_distance(HX, HZ) == 1

    A_terms = [(:x, 1), (:y, 1), (:y, 2)]
    B_terms = [(:x, 1), (:y, 1), (:y, 2)]
    HX_bb, HZ_bb = bb_check_matrices(2, 3, A_terms, B_terms)
    @test css_distance(HX_bb, HZ_bb) == 2

    @testset "simple hand-checkable example" begin
        # H =
        # [1 1 0 0
        #  1 0 1 0
        #  0 1 0 1]
        H = sparse(
            Int[1, 2, 1, 3, 2, 3],
            Int[1, 1, 2, 2, 3, 4],
            ones(Int, 6),
            3, 4,
        )

        check_to_bits, bit_to_checks, check_bit_pos, bit_check_pos = tanner_graph(H)

        # With SparseMatrixCSC / findnz, the traversal is column-major.
        @test check_to_bits == [
            Int[1, 2],
            Int[1, 3],
            Int[2, 4],
        ]

        @test bit_to_checks == [
            Int[1, 2],
            Int[1, 3],
            Int[2],
            Int[3],
        ]

        @test check_bit_pos[1] == Dict(1 => 1, 2 => 2)
        @test check_bit_pos[2] == Dict(1 => 1, 3 => 2)
        @test check_bit_pos[3] == Dict(2 => 1, 4 => 2)

        @test bit_check_pos[1] == Dict(1 => 1, 2 => 2)
        @test bit_check_pos[2] == Dict(1 => 1, 3 => 2)
        @test bit_check_pos[3] == Dict(2 => 1)
        @test bit_check_pos[4] == Dict(3 => 1)
    end

    @testset "isolated checks and bits" begin
        # H =
        # [0 1 0 0 0
        #  0 0 0 0 0
        #  0 0 0 1 0
        #  0 0 0 0 0]
        H = sparse(
            Int[1, 3],
            Int[2, 4],
            ones(Int, 2),
            4, 5,
        )

        check_to_bits, bit_to_checks, check_bit_pos, bit_check_pos = tanner_graph(H)

        @test check_to_bits == [
            Int[2],
            Int[],
            Int[4],
            Int[],
        ]

        @test bit_to_checks == [
            Int[],
            Int[1],
            Int[],
            Int[3],
            Int[],
        ]

        @test check_bit_pos[1] == Dict(2 => 1)
        @test isempty(check_bit_pos[2])
        @test check_bit_pos[3] == Dict(4 => 1)
        @test isempty(check_bit_pos[4])

        @test isempty(bit_check_pos[1])
        @test bit_check_pos[2] == Dict(1 => 1)
        @test isempty(bit_check_pos[3])
        @test bit_check_pos[4] == Dict(3 => 1)
        @test isempty(bit_check_pos[5])
    end

    @testset "generic consistency checks" begin
        # Intentionally not entered in sorted order.
        H = sparse(
            Int[3, 1, 2, 1, 3, 2],
            Int[2, 1, 3, 2, 4, 5],
            ones(Int, 6),
            4, 6,
        )

        check_to_bits, bit_to_checks, check_bit_pos, bit_check_pos = tanner_graph(H)
        M = Matrix(H)
        num_checks, num_bits = size(H)

        # Sizes match row/column weights.
        for i in 1:num_checks
            @test length(check_to_bits[i]) == count(!iszero, M[i, :])
            @test length(check_bit_pos[i]) == length(check_to_bits[i])
        end

        for j in 1:num_bits
            @test length(bit_to_checks[j]) == count(!iszero, M[:, j])
            @test length(bit_check_pos[j]) == length(bit_to_checks[j])
        end

        # Position dictionaries are exact inverses of adjacency lists.
        for i in 1:num_checks
            for (pos, j) in enumerate(check_to_bits[i])
                @test check_bit_pos[i][j] == pos
                @test M[i, j] != 0
            end
        end

        for j in 1:num_bits
            for (pos, i) in enumerate(bit_to_checks[j])
                @test bit_check_pos[j][i] == pos
                @test M[i, j] != 0
            end
        end

        # Every nonzero in H appears in both lookup tables, and no zero does.
        for i in 1:num_checks, j in 1:num_bits
            if M[i, j] != 0
                @test haskey(check_bit_pos[i], j)
                @test haskey(bit_check_pos[j], i)

                @test check_to_bits[i][check_bit_pos[i][j]] == j
                @test bit_to_checks[j][bit_check_pos[j][i]] == i
            else
                @test !haskey(check_bit_pos[i], j)
                @test !haskey(bit_check_pos[j], i)
            end
        end
    end    
end
