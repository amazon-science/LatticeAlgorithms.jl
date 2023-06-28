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
using BlockDiagonals

M = [1 -1 0 0;
     0 1 -1 0;
     0 0 1 -1;
     0 0 1  1]

lq_basis = lq_reduce(M)
L1, Q1 = lq_basis.L, lq_basis.Q

lll_basis = lll(M)
T2, L2, Q2 = lll_basis.T, lll_basis.L, lll_basis.Q

kz_basis = kz(M)
T3, L3, Q3 = kz_basis.T, kz_basis.L, kz_basis.Q


@test L1 * Q1 ≈ M
@test L1 ≈ [√2 0 0 0; -1/√2 √(3/2) 0 0; 0 -√(2/3) 2/√3 0; 0 -√(2/3) -1/√3 1]
@test Q1 ≈ [1/√2 -1/√2 0 0; 1/√6 1/√6 -√(2/3) 0; 1/(2√3) 1/(2√3) 1/(2√3) -√3/2; 1/2 1/2 1/2 1/2]
@test Q1 * transpose(Q1) ≈ Matrix(1I, size(Q1))
@test transpose(Q1) * Q1 ≈ Matrix(1I, size(Q1))
@test islowertriangular(L1)


@test T2 * L2 * Q2 ≈ M
@test islowertriangular(L2)
@test islllreduced(L2)
@test Q2 * transpose(Q2) ≈ Matrix(1I, size(Q2))
@test transpose(Q2) * Q2 ≈ Matrix(1I, size(Q2))
@test round.(Int, T2) ≈ T2
@test abs(det(T2)) ≈ 1

@test T3 * L3 * Q3 ≈ M
@test islowertriangular(L3)
# @test iskzreduced(L3)
@test Q3 * transpose(Q3) ≈ Matrix(1I, size(Q3))
@test transpose(Q3) * Q3 ≈ Matrix(1I, size(Q3))
@test round.(Int, T3) ≈ T3
@test abs(det(T3)) ≈ 1

@test_throws ErrorException lll(M, δ = 1.5)
@test_throws ErrorException kz(M, δ = 1.5)

# Massive test 


function test_reduction(N::Integer, method::String, samples::Integer=Int(1e2))
    elapsed_time = @elapsed for _ in 1 : samples

        M = rand(N, N) * N .- N/2

        if method == "LLL"
            basis = lll(M)
        elseif method == "KZ"
            basis = kz(M)
        elseif method == "LQ"
            basis = lq_reduce(M)
        end
        T, L, Q = basis.T, basis.L, basis.Q

        @test T * L * Q ≈ M
        @test islowertriangular(L)
        @test Q * transpose(Q) ≈ Matrix(1I, size(Q))
        @test transpose(Q) * Q ≈ Matrix(1I, size(Q))
        @test round.(Int, T) ≈ T
        @test abs(det(T)) ≈ 1
    end
    # println("test_reduction N=$N, elapsed_time=$elapsed_time")
end

for N in 1 : 6
    test_reduction(N, "LLL")
    test_reduction(N, "LQ")
    # test_reduction(N, "KZ")
end