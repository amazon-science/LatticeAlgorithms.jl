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

## Test cases for the old convention
# x1 = [-1.8, -1.4, 0.5, 0.0, 1, 2.1, 2.8, 3]
# f1 = [-2, -1, 0, 0, 1, 2, 3, 3] # closest points of x1
# g1 = [-1, -2, 1, 1, 2, 3, 2, 4] # second points of x1
# @test closest_integer.(x1) == f1
# @test second_closest_integer.(x1) == g1

## New convention for closest_integer and second_closest_integer
# Define the closest integers for some test numbers in [0, 1)
ks = [0.0, 0.1, 0.4, 0.5, 0.6, 0.9]
test_cases = Dict(
    0.0 => [0, 1], # The 1st and 2nd closest integers of 0.0
    0.9 => [1, 0], # The 1st and 2nd closest integers of 1.0
)
# The 1st and 2nd closest integers for numbers in [0.0 0.5] (inclusive for both ends) are the same as those for 0.0
test_cases[0.1] = test_cases[0.0] 
test_cases[0.4] = test_cases[0.0] 
test_cases[0.5] = test_cases[0.0] 

# The 1st and 2nd closest integers for numbers in (0.5, 1.0) (exclusive for both ends) are the same as those for 0.9
test_cases[0.6] = test_cases[0.9]

# Closest integers outside of [0, 1) are obtained by summing the integer part to the closest integers of its fractional part
# Define the closest integers for test numbers obtained by shifting the previous test points by some integers
for i in [-2, -1, 1, 2]
    for k in ks
        test_cases[k + i] = test_cases[k] .+ i
    end
end

# Run the tests
for (x, cps) in test_cases
    cps2 = [closest_integer(x), second_closest_integer(x)]
    @test cps == cps2
end

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