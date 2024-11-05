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

@test surface_code_Mp(3) == [
    1  0  0  1  0  0  0  0  0
    0  1  1  0  1  1  0  0  0
    0  0  2  0  0  0  0  0  0
    0  0  0  1  1  0  1  1  0
    0  0  0  0  2  0  0  0  0
    0  0  0  0  0  1  0  0  1
    0  0  0  0  0  0  2  0  0
    0  0  0  0  0  0  0  2  0
    0  0  0  0  0  0  0  0  2
]/√2

@test surface_code_Mq(3) == [
    1  1  0  1  1  0  0  0  0
    0  1  1  0  0  0  0  0  0
    0  0  2  0  0  0  0  0  0
    0  0  0  2  0  0  0  0  0
    0  0  0  0  1  1  0  1  1
    0  0  0  0  0  2  0  0  0
    0  0  0  0  0  0  1  1  0
    0  0  0  0  0  0  0  2  0
    0  0  0  0  0  0  0  0  2
]/√2

@test surface_code_M(3) == [
    1  0  1  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0
    0  1  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0
    0  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0
    0  0  0  1  0  1  0  0  0  1  0  1  0  0  0  0  0  0
    0  0  0  0  2  0  0  0  0  0  0  0  0  0  0  0  0  0
    0  0  0  0  0  2  0  0  0  0  0  0  0  0  0  0  0  0
    0  0  0  0  0  0  2  0  0  0  0  0  0  0  0  0  0  0
    0  0  0  0  0  0  0  1  0  1  0  0  0  1  0  1  0  0
    0  0  0  0  0  0  0  0  1  0  1  0  0  0  1  0  1  0
    0  0  0  0  0  0  0  0  0  2  0  0  0  0  0  0  0  0
    0  0  0  0  0  0  0  0  0  0  2  0  0  0  0  0  0  0
    0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  1
    0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0
    0  0  0  0  0  0  0  0  0  0  0  0  0  2  0  0  0  0
    0  0  0  0  0  0  0  0  0  0  0  0  0  0  2  0  0  0
    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2  0  0
    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2  0
    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2    
]/√2

@test abs(det(surface_code_M(3))) ≈ 2
@test abs(det(surface_code_M(5))) ≈ 2

@test surface_code_X_logicals(3)[1] == [1,4,7]
@test surface_code_Z_logicals(3)[1] == [1,2,3]

@test Set(values(surface_code_Z_stabilizers(3))) == Set([[2,3,5,6], [4,5,7,8], [1,4], [6,9]])
@test Set(values(surface_code_X_stabilizers(3))) == Set([[1,2,4,5], [5,6,8,9], [2,3], [7,8]])
@test Set(values(surface_code_stabilizers(3))) == Set([[2, 8],[9, 11, 15, 17], [4, 6, 10, 12], [8, 10, 14, 16], [13, 15], [12, 18], [3, 5], [1, 3, 7, 9]])

@test Set(values(surface_code_Z_stabilizers(5))) == Set(
    [
        [2,3,7,8], [4,5,9,10], [6,7,11,12], [8,9,13,14], [1,6], [11, 16],
        [12,13,17,18], [14,15,19,20], [16,17,21,22], [18,19,23,24], [10, 15], [20, 25]
    ]
)

@test Set(values(surface_code_X_stabilizers(5))) == Set(
    [
        [1,2,6,7], [3,4,8,9], [7,8,12,13], [9,10,14,15], [2,3], [4,5],
        [11,12,16,17], [13,14,18,19], [17,18,22,23], [19,20,24,25], [21,22], [23,24]
    ]
)

function test_surface_code_stabilizers(d)
    
    # Test surface_code_stabilizers via generator matrix
    M1 = surface_code_M(d)

    M2 = 2 * diagm(ones(Int64, 2d^2))
    dict = surface_code_stabilizers(d)
    for (_, value_list) in dict
        min_value = minimum(value_list)
        for value in value_list
            M2[min_value, value] = 1
        end
    end
    M2 = M2/√2

    @test M1 ≈ M2
end

for d in 3:2:19
    test_surface_code_stabilizers(d)
end


drange = [3,5]
num_samples = 1e3
for d in drange
    x_stab_generators = collect(values(surface_code_X_stabilizers(d)))
    x_stabilizers = get_stabilizer_group_from_generators(collect(values(x_stab_generators)))
        
    z_stab_generators = collect(values(surface_code_Z_stabilizers(d)))
    z_stabilizers = get_stabilizer_group_from_generators(collect(values(z_stab_generators)))

    for (whichtype, stabilizers) in zip(["x", "z"], [x_stabilizers, z_stabilizers])
        println("d, whichtype =$d, $whichtype")
        
        for _ in 1 : num_samples
            ws = 0.6 * randn(d^2)
            ws = abs.(ws)
            Zw = sum([prod(ws[stab]) for stab in stabilizers])
            
            ws_unrotated = LatticeAlgorithms.embed_weights(ws, whichtype)
            Zw2 = bsv(ws_unrotated)

            @assert Zw2 ≈ log10(Zw)
        end
    end
end


num_samples = 1e3
σ = 0.6

for d in [3, 5]
    num_qubits = d^2
    stab = surface_code_X_stabilizers(d)
    stabilizers = get_stabilizer_group_from_generators(collect(values(stab)))

    for _ in 1 : num_samples
        ηs = σ * randn(num_qubits)
        ws = get_ws_square(ηs, σ)
        w2s = get_ws_square(ηs .+ √π, σ)
        
        expected_value = 0
        for stab in stabilizers
            prodd = [i ∈ stab ? ws[i] : w2s[i] for i in 1 : length(ηs)]
            expected_value += prod(prodd)
        end
        
        w3s = ws ./ w2s
        w3s_unrotated = LatticeAlgorithms.embed_weights(w3s, "x")
        value = bsv(w3s_unrotated) + log10(prod(w2s))
        
        if !(10^value ≈ expected_value)
            println(ws)
            println(10^value)
            println(expected_value)
            println()
        end
        @assert value ≈ log10(expected_value)
    end
end


num_samples = 1e3
σ = 0.6

for d in [3, 5]
    num_qubits = d^2
    
    stab_x = surface_code_X_stabilizers(d)
    stabilizers_x = get_stabilizer_group_from_generators(collect(values(stab_x)))
    x_logical = surface_code_X_logicals(d)[1]

    stab_z = surface_code_Z_stabilizers(d)
    stabilizers_z = get_stabilizer_group_from_generators(collect(values(stab_z)))
    z_logical = surface_code_Z_logicals(d)[1]

    for _ in 1 : num_samples
        ηs = σ * randn(num_qubits)

        # test for x 
        rec = bsv_surface_code(ηs, σ; Nv = 5, subspace="x")
        lstar2, _, _ = brute_force_mld_concatenated_square(ηs, σ, x_logical, stabilizers_x)
        rec2 = -(ηs - lstar2)
        @test rec ≈ rec2

        # test for z
        rec = bsv_surface_code(ηs, σ; Nv = 5, subspace="z")
        lstar2, _, _ = brute_force_mld_concatenated_square(ηs, σ, z_logical, stabilizers_z)
        rec2 = -(ηs - lstar2)
        @test rec ≈ rec2

    end
end

# Test for tn_surf_hex

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
    num_qubits = d^2
    bigX = vcat([X for _ in 1 : num_qubits]...)
    bigZ = vcat([Z for _ in 1 : num_qubits]...)
    stabs = surface_code_stabilizers(d)
    full_stabs = get_stabilizer_group_from_generators(collect(values(stabs)))

    indicators = get_indicators_from_stabilizers(full_stabs)
    TN, indices = tn_template_surf_hex(d)
    for _ in 1 : num_samples
        ηs = σ * randn(2num_qubits)
        lstar, prob_I, prob_X, prob_Y, prob_Z, _, _, _, _ = tn_surf_hex(ηs, σ, TN, indices, bigZ, bigX, χ; S=S, Nv=Nv)

        lstar2, prob_I2, prob_X2, prob_Y2, prob_Z2 = brute_force_mld_concatenated_non_square(ηs, σ, bigX, bigZ, indicators; S=S, Nv=Nv)

        @test lstar == lstar2
        @test log10(prob_I2) ≈ prob_I
        @test log10(prob_X2) ≈ prob_X
        @test log10(prob_Y2) ≈ prob_Y
        @test log10(prob_Z2) ≈ prob_Z
    end
end
