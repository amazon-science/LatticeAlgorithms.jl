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
using LinearAlgebra
using Test

## Triangular color codes

@test triangular_color_code_num_qubits(3) == 7
@test triangular_color_code_num_qubits(5) == 19

@test Set(values(triangular_color_code_stabilizers(3, "R"))) == Set([[2,3,5,6]])
@test Set(values(triangular_color_code_stabilizers(3, "G"))) == Set([[1,2,3,4]])
@test Set(values(triangular_color_code_stabilizers(3, "B"))) == Set([[3,4,6,7]])

@test Set(values(triangular_color_code_stabilizers(5, "R"))) == Set([[2,3,5,6], [8, 11, 15, 16], [9, 10, 12, 13, 17, 18]])
@test Set(values(triangular_color_code_stabilizers(5, "G"))) == Set([[1,2,3,4], [5, 6, 8, 9, 11, 12], [7, 10, 13, 14]])
@test Set(values(triangular_color_code_stabilizers(5, "B"))) == Set([[3,4,6,7,9, 10], [11, 12, 16, 17], [13, 14, 18, 19]])

@test Set(values(triangular_color_code_stabilizers(3))) == Set([[2,3,5,6], [1,2,3,4], [3,4,6,7]])

@test Set(values(triangular_color_code_stabilizers(5))) == Set([[3,4,6,7,9, 10], [11, 12, 16, 17], [13, 14, 18, 19], [1,2,3,4], [5, 6, 8, 9, 11, 12], [7, 10, 13, 14], [2,3,5,6], [8, 11, 15, 16], [9, 10, 12, 13, 17, 18]])

@test Set(values(triangular_color_code_logicals(3))) == Set([[1, 2, 5]])
@test Set(values(triangular_color_code_logicals(5))) == Set([[1, 2, 5, 8, 15]])

@test Set(values(triangular_color_code_normalizers(3))) == Set([[1, 2, 5], [2,3,5,6], [1,2,3,4], [3,4,6,7]])

@test Set(values(triangular_color_code_normalizers(5))) == Set([[1, 2, 5, 8, 15], [3,4,6,7,9, 10], [11, 12, 16, 17], [13, 14, 18, 19], [1,2,3,4], [5, 6, 8, 9, 11, 12], [7, 10, 13, 14], [2,3,5,6], [8, 11, 15, 16], [9, 10, 12, 13, 17, 18]])

@test triangular_color_code_Mq(3) == [
    1  1  1  1  0  0  0
    0  1  1  0  1  1  0
    0  0  1  1  0  1  1
    0  0  0  2  0  0  0
    0  0  0  0  2  0  0
    0  0  0  0  0  2  0
    0  0  0  0  0  0  2
]/√2

@test triangular_color_code_Mq(3) == triangular_color_code_Mp(3)

@test abs(det(triangular_color_code_M(3))) ≈ 2
@test abs(det(triangular_color_code_M(5))) ≈ 2

# function test_triangular_color_code_distances(d)
#     ds = distances(triangular_color_code_M(d))
#     @test ds[1]≈ds[2]
#     @test ds[1]≈ds[3]
# end

# for d in [3]
#     test_triangular_color_code_distances(d)
# end

## Octagonal color codes

@test octagonal_color_code_num_qubits(2) == 8
@test octagonal_color_code_num_qubits(4) == 32
@test octagonal_color_code_num_qubits(6) == 72


@test Set(values(octagonal_color_code_stabilizers(2, "R"))) == Set([[1,2,7,8], [3,4,5,6]])
@test Set(values(octagonal_color_code_stabilizers(2, "G"))) == Set([[1, 2, 3, 4, 7, 8, 5, 6]])
@test Set(values(octagonal_color_code_stabilizers(2, "B"))) == Set([[1, 3, 5, 7, 2, 4, 6, 8]])

@test Set(values(octagonal_color_code_stabilizers(4, "R"))) == Set([[21, 22, 25, 26],
[7, 8, 11, 12],
[5, 6, 9, 10],
[23, 24, 27, 28],
[2, 3, 30, 31],
[14, 15, 18, 19],
[13, 16, 17, 20],
[1, 4, 29, 32]])

@test Set(values(octagonal_color_code_stabilizers(4, "G"))) == Set([[9, 10, 13, 14, 17, 18, 21, 22], [1, 2, 5, 6, 29, 30, 25, 26], [3, 4, 7, 8, 31, 32, 27, 28], [11, 12, 15, 16, 19, 20, 23, 24]])

@test Set(values(octagonal_color_code_stabilizers(4, "B"))) == Set([[17, 21, 25, 29, 20, 24, 28, 32], [2, 3, 6, 7, 10, 11, 14, 15], [1, 5, 9, 13, 4, 8, 12, 16], [18, 19, 22, 23, 26, 27, 30, 31]])

@test Set(values(octagonal_color_code_logicals(2))) == Set([[1, 2], [3, 4], [3, 5], [1, 7]])
@test Set(values(octagonal_color_code_logicals(4))) == Set([[1,2,3,4], [5,6,7,8], [1,13,17,29], [5,9,21,25]])

@test Set(values(octagonal_color_code_normalizers(2))) == Set([[1, 2], [3, 4], [3, 5], [1, 7], [1,2,7,8], [3,4,5,6]])
# @test Set(values(octagonal_color_code_normalizers(4))) == Set([[1,2,3,4], [5,6,7,8], [1,13,17,29], [5,9,21,25]])

@test octagonal_color_code_Mp(2) == [
    1  1  0  0  0  0  1  1
    0  2  0  0  0  0  0  0
    0  0  1  1  1  1  0  0
    0  0  0  2  0  0  0  0
    0  0  0  0  2  0  0  0
    0  0  0  0  0  2  0  0
    0  0  0  0  0  0  2  0
    0  0  0  0  0  0  0  2
]/√2

@test octagonal_color_code_Mp(4) == octagonal_color_code_Mp(4)
@test octagonal_color_code_Mp(2) == octagonal_color_code_Mp(2)

# Test for tn_color_square

num_samples = 1e3
σ = 0.6
Nv = 5
χ = 100

for d in [3, 5]
    num_qubits = triangular_color_code_num_qubits(d)
    stabs = triangular_color_code_stabilizers(d)
    full_stabs = get_stabilizer_group_from_generators(collect(values(stabs)))

    indicators = get_indicators_from_stabilizers(full_stabs)
    TN, indices = tn_template_color_square(d)

    logical = triangular_color_code_logicals(d)[1]
    for _ in 1 : num_samples
        ηs_q = σ * randn(num_qubits)
        lstar, prob_I, prob_X, _, _ = tn_color_square(ηs_q, σ, logical, TN, indices, χ; Nv=Nv)

        lstar2, prob_I2, prob_X2 = brute_force_mld_concatenated_square(ηs_q, σ, logical, full_stabs; Nv=Nv)

        @test lstar == lstar2
        @test log10(prob_I2) ≈ prob_I
        @test log10(prob_X2) ≈ prob_X
    end
end


# Test for tn_color_hex

num_samples = 1e3
σ = 0.6
Nv = 5
χ = 100
S = [2 1; 0 sqrt(3)] / (12)^(1/4)
S_T = transpose(S)
X = √π * S * [1, 0]
Z = √π * S * [0, 1]
Y = X + Z

for d in [3]
    num_qubits = triangular_color_code_num_qubits(d)
    bigX = vcat([X for _ in 1 : num_qubits]...)
    bigZ = vcat([Z for _ in 1 : num_qubits]...)
    
    stabs0 = triangular_color_code_stabilizers(d)
    stabs1 = Dict{Int64, Vector{Int64}}()
    for (k, v) in stabs0
        stabs1[k] = 2 .* v .- 1
        stabs0[k] = 2 .* v
    end
    stabs = merge_two_stabilizers(stabs0, stabs1)

    full_stabs = get_stabilizer_group_from_generators(collect(values(stabs)))

    indicators = get_indicators_from_stabilizers(full_stabs)
    TN, indices = tn_template_color_hex(d)
    for _ in 1 : num_samples
        ηs = σ * randn(2num_qubits)
        lstar, prob_I, prob_X, prob_Y, prob_Z, _, _, _, _ = tn_color_hex(ηs, σ, TN, indices, bigZ, bigX, χ; S=S, Nv=Nv)

        lstar2, prob_I2, prob_X2, prob_Y2, prob_Z2 = brute_force_mld_concatenated_non_square(ηs, σ, bigX, bigZ, indicators; S=S, Nv=Nv)

        @test lstar == lstar2
        @test log10(prob_I2) ≈ prob_I
        @test log10(prob_X2) ≈ prob_X
        @test log10(prob_Y2) ≈ prob_Y
        @test log10(prob_Z2) ≈ prob_Z
    end
end
