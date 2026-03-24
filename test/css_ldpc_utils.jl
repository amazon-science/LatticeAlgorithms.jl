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
end
