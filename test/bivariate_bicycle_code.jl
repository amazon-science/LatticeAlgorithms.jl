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

function test_bb_code_basic_properties(l, m, A_terms, B_terms, n_expected, k_expected)
    @test bb_num_data_qubits(l, m) == n_expected
    @test bb_num_logicals(l, m, A_terms, B_terms) == k_expected
    @test bb_parameters(l, m, A_terms, B_terms) == (n=n_expected, k=k_expected)

    A, B = bb_matrices(l, m, A_terms, B_terms)
    HX, HZ = bb_check_matrices(l, m, A_terms, B_terms)

    @test size(A) == (l * m, l * m)
    @test size(B) == (l * m, l * m)
    @test size(HX) == (l * m, 2 * l * m)
    @test size(HZ) == (l * m, 2 * l * m)

    @test HX == hcat(A, B)
    @test HZ == hcat(transpose(B), transpose(A))

    # CSS commutation condition
    @test all(x -> x == 0, mod.(HX * transpose(HZ), 2))

    # All Table-3 BB examples have weight-6 checks.
    x_stabs = bb_X_stabilizers(l, m, A_terms, B_terms)
    z_stabs = bb_Z_stabilizers(l, m, A_terms, B_terms)

    @test length(x_stabs) == l * m
    @test length(z_stabs) == l * m
    @test all(length(v) == 6 for v in values(x_stabs))
    @test all(length(v) == 6 for v in values(z_stabs))

    # Combined stabilizer dictionary convention:
    # X on qubit j -> 2j-1, Z on qubit j -> 2j
    stabs = bb_stabilizers(l, m, A_terms, B_terms)
    @test length(stabs) == 2 * l * m

    for i in 1:(l * m)
        @test haskey(stabs, i)
        @test all(isodd, stabs[i])
    end

    for i in (l * m + 1):(2 * l * m)
        @test haskey(stabs, i)
        @test all(iseven, stabs[i])
    end

    # Logical operators
    X_logs = bb_X_logicals(l, m, A_terms, B_terms)
    Z_logs = bb_Z_logicals(l, m, A_terms, B_terms)
    @test length(X_logs) == k_expected
    @test length(Z_logs) == k_expected

    # GKP generators
    Mq = bb_Mq(l, m, A_terms, B_terms)
    Mp = bb_Mp(l, m, A_terms, B_terms)
    M = bb_M(l, m, A_terms, B_terms)

    @test size(Mq) == (n_expected, n_expected)
    @test size(Mp) == (n_expected, n_expected)
    @test size(M) == (2 * n_expected, 2 * n_expected)
end

# Table-3 examples from the BB-code paper
examples = [
    (
        6, 6,
        [(:x, 3), (:y, 1), (:y, 2)],
        [(:y, 3), (:x, 1), (:x, 2)],
        72, 12,
    ),
    (
        15, 3,
        [(:x, 9), (:y, 1), (:y, 2)],
        [(:I, 0), (:x, 2), (:x, 7)],
        90, 8,
    ),
    (
        9, 6,
        [(:x, 3), (:y, 1), (:y, 2)],
        [(:y, 3), (:x, 1), (:x, 2)],
        108, 8,
    ),
    (
        12, 6,
        [(:x, 3), (:y, 1), (:y, 2)],
        [(:y, 3), (:x, 1), (:x, 2)],
        144, 12,
    ),
    (
        12, 12,
        [(:x, 3), (:y, 2), (:y, 7)],
        [(:y, 3), (:x, 1), (:x, 2)],
        288, 12,
    ),
]

for (l, m, A_terms, B_terms, n_expected, k_expected) in examples
    test_bb_code_basic_properties(l, m, A_terms, B_terms, n_expected, k_expected)
end

# Basic sanity checks for exponent reduction
let l = 6, m = 6
    x = bb_x_matrix(l, m)
    y = bb_y_matrix(l, m)
    @test x^l == Matrix{Int64}(I, l * m, l * m)
    @test y^m == Matrix{Int64}(I, l * m, l * m)
end

# Input validation
@test_throws ErrorException bb_x_matrix(0, 3)
@test_throws ErrorException bb_y_matrix(3, 0)
@test_throws ErrorException LatticeAlgorithms._bb_matrix(6, 6, [(:x, 1), (:y, 1)])
@test_throws ErrorException LatticeAlgorithms._bb_matrix(6, 6, [(:x, 1), (:x, 1), (:y, 2)])

HX, HZ = bb_check_matrices(6, 6,
    [(:x, 3), (:y, 1), (:y, 2)],
    [(:y, 3), (:x, 1), (:x, 2)]
)

@test css_distance(HX, HZ; sector=:Z) == 6
@test css_distance(HX, HZ; sector=:X) == 6