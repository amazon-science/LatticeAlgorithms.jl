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

σ = 0.6
Nv = 1
S = [2 1; 0 sqrt(3)] / (12)^(1/4)
S_T = transpose(S)
X = √π * S * [1, 0]
Z = √π * S * [0, 1]
Y = X + Z


ws = get_ws_square([0.1], σ; Nv=Nv)
@test ws ≈ [exp(-(0.1-√π)^2/(2 * 0.6^2))]

x = [0.1, 0.1]
ws = get_ws_non_square(x, σ; Nv=Nv)
cp_I = closest_point(x, 2√π * S_T)
cp_X = closest_point(x-X, 2√π * S_T)
cp_Y = closest_point(x-Y, 2√π * S_T)
cp_Z = closest_point(x-Z, 2√π * S_T)
@test ws ≈ [
    [
        exp(-norm(x-cp_I)^2/(2σ^2)),
        exp(-norm(x-X-cp_X)^2/(2σ^2)),
        exp(-norm(x-Z-cp_Z)^2/(2σ^2)),
        exp(-norm(x-Y-cp_Y)^2/(2σ^2))
    ]
]



d = 3
logical = triangular_color_code_logicals(d)[1]
stabilizers = triangular_color_code_stabilizers(d)
full_stabilizers = get_stabilizer_group_from_generators(collect(values(stabilizers)))
num_qubits = triangular_color_code_num_qubits(d)
lstar, _ , _ = brute_force_mld_concatenated_square([0 for _ in 1 : num_qubits], σ, logical, full_stabilizers)

@test lstar ≈ [0.0 for _ in 1 : num_qubits]

bigS_T = BlockDiagonal([S_T for _ in 1 : num_qubits])

bigX = vcat([X for _ in 1 : num_qubits]...)
bigZ = vcat([Z for _ in 1 : num_qubits]...)

indicators = get_indicators_from_stabilizers(full_stabilizers)
@test Set(indicators) == Set([
    [0, 0, 0, 0, 0, 0, 0],
    [2, 1, 3, 0, 0, 0, 0],
    [0, 3, 2, 1, 0, 0, 0],
    [2, 2, 1, 1, 0, 0, 0],
    [3, 3, 0, 0, 0, 0, 0],
    [1, 2, 3, 0, 0, 0, 0],
    [3, 0, 2, 1, 0, 0, 0],
    [1, 1, 1, 1, 0, 0, 0],   
])

lstar, _, _, _, _ = brute_force_mld_concatenated_non_square([0.0 for _ in 1 : 2num_qubits], σ, bigX, bigZ, indicators; S=S)

@test lstar ≈ [0.0 for _ in 1 : 2num_qubits]