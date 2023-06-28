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

function test_decode_rep_rec(N, num_samples=1e2)
    M = rep_rec(N)
    Mperp = GKP_logical_operator_generator_canonical(M)
    elapsed_time = @elapsed for _ in 1 : num_samples
        xs = [rand(2N) * 2N .- 2N/2 for _ in 1 : num_samples]
        y1s = [closest_point(x, √(2π) * Mperp) for x in xs]
        y2s = [decode_rep_rec(x) for x in xs]
        @test y1s ≈ y2s
    end
    # println("test_decode_rep_rec for N = $N, elapsed_time=$elapsed_time")
end

for N in 2 : 6
    test_decode_rep_rec(N)
end

function test_decode_YY_rep_rec(N, num_samples=1e2)
    M = YY_rep_rec(N)
    Mperp = GKP_logical_operator_generator_canonical(M)
    _, _, Q2, r2 = tlq_YY_rep_rec(N)

    elapsed_time = @elapsed for _ in 1 : num_samples
        xs = [rand(4N) * 4N .- 4N/2 for _ in 1 : num_samples]
        y1s = [closest_point(x, √(2π) * Mperp) for x in xs]
        y2s = [decode_YY_rep_rec(x, Q2, r2) for x in xs]
        @test y1s ≈ y2s
    end
    # println("test_decode_YY_rep_rec for N = $N, elapsed_time=$elapsed_time")
end

for N in 2 : 6
    test_decode_YY_rep_rec(N)
end

