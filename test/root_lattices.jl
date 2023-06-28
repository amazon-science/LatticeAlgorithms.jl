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

x1 = [-1.8, -1.4, 0.5, 0.0, 1, 2.1, 2.8, 3]
f1 = [-2, -1, 0, 0, 1, 2, 3, 3] # closest points of x1
g1 = [-1, -2, 1, 1, 2, 3, 2, 4] # second points of x1

@test closest_integer.(x1) == f1

@test second_closest_integer.(x1) == g1

@test Dn(2) == [1 -1; 1 1]
@test Dn(3) == [1 -1 0; 0 1 -1; 0 1 1]

@test Dn_dual(2) == Matrix(inv(transpose([1 -1; 1 1])))
@test Dn_dual(3) == Matrix(inv(transpose([1 -1 0; 0 1 -1; 0 1 1])))

## massive tests

function test_closest_point_Zn(N, num_samples=1e2)
    elapsed_time = @elapsed for _ in 1 : num_samples
        xs = [rand(N) * N .- N/2 for _ in 1 : num_samples]
        us = [closest_point(x, Matrix(1I, N,N)) for x in xs]
        fxs = [closest_point_Zn(x) for x in xs]
        @test us ≈ fxs
    end
    # println("test_closest_point_Zn for N = $N, elapsed_time=$elapsed_time")
end

for N in 2 : 10
    test_closest_point_Zn(N)
end


function test_closest_point_scaled_Zn_1(N, num_samples=1e2)
    elapsed_time = @elapsed for _ in 1 : num_samples
        scale = rand() * N - N/2
        Zn_scaled = scale * Matrix(1I, N, N)
        xs = [rand(N) * N .- N/2 for _ in 1 : num_samples]
        us = [closest_point(x, Zn_scaled) for x in xs]
        fxs = [closest_point_scaled_Zn(x, scale) for x in xs]
        @test us ≈ fxs
    end
    # println("test_closest_point_scaled_Zn_1 for N = $N, elapsed_time=$elapsed_time")
end

for N in 2 : 10
    test_closest_point_scaled_Zn_1(N)
end


function test_closest_point_scaled_Zn_2(N, num_samples=1e2)
    elapsed_time = @elapsed for _ in 1 : num_samples
        scales = rand(N) * N .- N/2
        Zn_scaled = diagm(vec(scales))
        xs = [rand(N) * N .- N/2 for _ in 1 : num_samples]
        us = [closest_point(x, Zn_scaled) for x in xs]
        fxs = [closest_point_scaled_Zn(x, scales) for x in xs]
        @test us ≈ fxs
    end
    # println("test_closest_point_scaled_Zn_2 for N = $N, elapsed_time=$elapsed_time")
end

for N in 2 : 10
    test_closest_point_scaled_Zn_2(N)
end


function test_closest_point_scaled_Dn_1(N, num_samples=1e2)
    elapsed_time = @elapsed for _ in 1 : num_samples
        scales = rand(N) * N .- N/2
        DN = Dn(N)
        DN_scaled = DN * diagm(vec(scales))
        xs = [rand(N) * N .- N/2 for _ in 1 : num_samples]
        u0 = [closest_point_scaled_Dn(x, scales) for x in xs]
        u2 = [closest_point(x, DN_scaled) for x in xs]

        @test u0 ≈ u2 
    end
    # println("test_closest_point_scaled_Dn_1 for N = $N, elapsed_time=$elapsed_time")
end

for N in 2 : 10
    test_closest_point_scaled_Dn_1(N)
end


function test_closest_point_scaled_Dn_2(N, num_samples=1e2)
    elapsed_time = @elapsed for _ in 1 : num_samples
        scale = rand() * N - N/2
        DN = Dn(N)
        DN_scaled = DN * scale
        xs = [rand(N) * N .- N/2 for _ in 1 : num_samples]
        u0 = [closest_point_scaled_Dn(x, scale) for x in xs]
        u2 = [closest_point(x, DN_scaled) for x in xs]

        @test u0 ≈ u2 
    end
    # println("test_closest_point_scaled_Dn_2 for N = $N, elapsed_time=$elapsed_time")
end

for N in 2 : 10
    test_closest_point_scaled_Dn_2(N)
end

function test_closest_point_Dn(N, num_samples=1e2)
    elapsed_time = @elapsed for _ in 1 : num_samples
        DN = Dn(N)
        xs = [rand(N) * N .- N/2 for _ in 1 : num_samples]
        u0 = [closest_point_Dn(x) for x in xs]
        u2 = [closest_point(x, DN) for x in xs]

        @test u0 ≈ u2 
    end
    # println("test_closest_point_Dn for N = $N, elapsed_time=$elapsed_time")
end

for N in 2 : 10
    test_closest_point_Dn(N)
end



function test_closest_point_scaled_Dn_dual_1(N, num_samples=1e2)
    elapsed_time = @elapsed for _ in 1 : num_samples
        scales = rand(N) * N .- N/2
        DN_dual = Dn_dual(N)
        DN_dual_scaled = DN_dual * diagm(vec(scales))
        xs = [rand(N) * N .- N/2 for _ in 1 : num_samples]
        u0 = [closest_point_scaled_Dn_dual(x, scales) for x in xs]
        u2 = [closest_point(x, DN_dual_scaled) for x in xs]

        @test u0 ≈ u2 
    end
    # println("test_closest_point_scaled_Dn_dual_1 for N = $N, elapsed_time=$elapsed_time")
end

for N in 2 : 10
    test_closest_point_scaled_Dn_dual_1(N)
end

function test_closest_point_scaled_Dn_dual_2(N, num_samples=1e2)
    elapsed_time = @elapsed for _ in 1 : num_samples
        scale = rand() * N - N/2
        DN_dual = Dn_dual(N)
        DN_dual_scaled = DN_dual * scale
        xs = [rand(N) * N .- N/2 for _ in 1 : num_samples]
        u0 = [closest_point_scaled_Dn_dual(x, scale) for x in xs]
        u2 = [closest_point(x, DN_dual_scaled) for x in xs]

        @test u0 ≈ u2 
    end
    # println("test_closest_point_scaled_Dn_dual_2 for N = $N, elapsed_time=$elapsed_time")
end

for N in 2 : 10
    test_closest_point_scaled_Dn_dual_2(N)
end

function test_closest_point_Dn_dual(N, num_samples=1e2)
    elapsed_time = @elapsed for _ in 1 : num_samples
        DN_dual = Dn_dual(N)
        xs = [rand(N) * N .- N/2 for _ in 1 : num_samples]
        u0 = [closest_point_Dn_dual(x) for x in xs]
        u2 = [closest_point(x, DN_dual) for x in xs]

        @test u0 ≈ u2 
    end
    # println("test_closest_point_Dn_dual for N = $N, elapsed_time=$elapsed_time")
end

for N in 2 : 10
    test_closest_point_Dn_dual(N)
end



function test_An_and_An_dual(N)
    @test An(N) * transpose(An(N)) ≈ Matrix(1I, N, N) + ones(N, N)
    @test An_dual(N) * transpose(An_dual(N)) ≈ Matrix(1I, N, N) - ones(N, N)/(N+1)
    # println("test_An for N = $N, elapsed_time=$elapsed_time")
end

for N in 2 : 10
    test_An_and_An_dual(N)
end

@test det(E8()) ≈ 1