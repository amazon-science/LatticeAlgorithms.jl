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


function test_decode_concatenated_rectangular_GKP_code(N, num_samples=1e4)

    η = rand(N)
    ηs = vcat([[ηη, 1/ηη] for ηη in η]...)
    
    M = 2 * Matrix(1I, 2N, 2N)    
    stab_rep_rec = Dict{Int64, Vector{Int64}}()
    for i = 1 : N-1
        M[2i-1, 2i-1:2i+2] = [1 0 1 0]
        stab_rep_rec[i] = [i, i+1]
    end
    M = M/√2 * diagm(ηs)
    
    Mperp_rep_rec = GKP_logical_operator_generator_canonical(M)
    
    elapsed_time = @elapsed for _ in 1 : num_samples
        ξ = 10 * rand(2N)
        y = closest_point(ξ, √(2π) * Mperp_rep_rec)
        y2 = decode_concatenated_rectangular_GKP_code(ξ[2:2:end], stab_rep_rec, ηs[2:2:end])
    
        @test max((abs.(y[2:2:end] - y2[1]))...) < 1e-10
    end

    # println("test_decode_concatenated_rectangular_GKP_code for d = $d, elapsed_time=$elapsed_time")
end

for N in [3,5]
    test_decode_concatenated_rectangular_GKP_code(N)
end



function decode_concatenated_rectangular_GKP_code_non_CSS_surf(d, num_samples=1e2)
    M = surface_code_M(d)
    Mperp = GKP_logical_operator_generator(M)

    surface_code_x_stabilizers = surface_code_X_stabilizers(d)
    surface_code_z_stabilizers = surface_code_Z_stabilizers(d)
    surf_stabilizers = surface_code_stabilizers(d)

    elapsed_time = @elapsed for _ in 1 : num_samples
        xs = [rand(2d^2) * 2d^2 .- 2d^2/2 for _ in 1 : num_samples]
        y1s = [closest_point(x, √(2π) * Mperp) for x in xs]
        y2s = [decode_concatenated_rectangular_GKP_code(x, surf_stabilizers, ones(length(x)), type_stabilizers="non-CSS")[1] for x in xs]
        @test y1s ≈ y2s
    end
    # println("decode_concatenated_rectangular_GKP_code_non_CSS_surf for d = $d, elapsed_time=$elapsed_time")
end

for d in [3]
    decode_concatenated_rectangular_GKP_code_non_CSS_surf(d)
end

function decode_concatenated_rectangular_GKP_code_non_CSS_YY_rep_rec(N, num_samples=1e2)
    xs = [rand(4N) * 4N .- 4N/2 for _ in 1 : num_samples]

    # First way to get the generator matrix for the YY_rep_rec code
    M0 = YY_rep_rec(N)
    Mperp0 = GKP_logical_operator_generator_canonical(M0)
    y1s = [closest_point(x, √(2π) * Mperp0) for x in xs]
    
    # Second way to get the generator matrix for the YY_rep_rec code
    # for confirming that the stabilizers are correct.
    YY_rep_rec_stabilizers = Dict{Int64, Vector{Int64}}()
    M = 2 * diagm(ones(Int64, 4N))
    for i in 1 : N-1
        YY_rep_rec_stabilizers[i] = [2i-1, 2i+1]
        M[2i-1, 2i-1] = M[2i-1, 2i+1] = 1
    end
    for i in 1 : N-1
        YY_rep_rec_stabilizers[N-1+i] = [2N + 2i-1, 2N + 2i+1]
        M[2N + 2i-1, 2N + 2i-1] = M[2N + 2i-1, 2N + 2i+1] = 1
    end
    YY_rep_rec_stabilizers[2N-1] = sort(vcat([1, 2N+1], 2*(1:2N)))
    M[end, 1] = M[end, 2N+1] = 1
    M[end, 2*(1:2N)] .= 1
    Z = vcat([[N^(1/4), N^(-1/4)] for _ in 1 : 2N]...)
    M = M/√2 * diagm(Z)
    
    Mperp = GKP_logical_operator_generator_canonical(M)
    y2s = [closest_point(x, √(2π) * Mperp) for x in xs]

    # Third way to get the closest point for the YY_rep_rec code
    # using the non-CSS stabilizers
    y3s = [decode_concatenated_rectangular_GKP_code(x, 
            YY_rep_rec_stabilizers, 
            Z, 
            type_stabilizers="non-CSS")[1] 
           for x in xs
    ]

    @test y1s≈y2s
    @test y1s≈y3s
end

for N in [2,3,4,5]
    decode_concatenated_rectangular_GKP_code_non_CSS_YY_rep_rec(N)
end

