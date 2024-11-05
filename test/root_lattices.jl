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


####### Below for K closest points

@test 1 == heaviside(0.1) == heaviside(1.0)
@test 0 == heaviside(0) == heaviside(0.0) == heaviside(-0.0) == heaviside(-0.1)

### Test next_closest_integer
# Define the closest points for some test points in [0, 1)
ks = [0.0, 0.1, 0.4, 0.5, 0.6, 0.9]
test_cases = Dict(
    0.0 => [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5], # The first few closest integers of 0.0
    0.9 => [1, 0, 2, -1, 3, -2, 4, -3, 5, -4, 6, -5], # The first few closest integers of 0.9
)

# The closest integers for numbers in [0.0 0.5] (inclusive for both ends) are the same as those for 0.0
test_cases[0.1] = test_cases[0.0]
test_cases[0.4] = test_cases[0.0]
test_cases[0.5] = test_cases[0.0]
# The closest integers for numbers in (0.5, 1.0) (exclusive for both ends) are the same as those for 0.9
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
    cp = cps[1]
    cps2 = [cp]
    for _ in 2 : length(cps)
        w = next_closest_integer(x, cps2[end])
        push!(cps2, w)
    end
    @test cps == cps2
end

## Define some test cases for closest_points_Zn

# For the 1D Zn lattice, the results should be the same as closest_integers (or next_closest_integer)
for (x, cps) in test_cases
    # println(x)
    cps2 = closest_points_Zn([x], length(cps))
    cps2 = [cp[1] for cp in cps2]
    @test cps == cps2
end

## 2D and 3D for random x


"""
    all_closest_points_within_distance(x::Vector, dist::Real)

Return the all the closest points in Zn to x within a given distance dist.

Note: This is for testing purpose only and only dimension 2 and 3 are implemented.
"""
function all_closest_points_within_distance(x::Vector, dist::Real)

    if length(x) ∉ [2, 3]
        error("only implemented for dimension 2 and 3 for testing.")
    end

    cps = []
    dists = []
    if length(x) == 2
        x2s = Int.(floor(x[1] - dist) : 1 : ceil(x[1] + dist))
        y2s = Int.(floor(x[2] - dist) : 1 : ceil(x[2] + dist))

        for x2 in x2s
            for y2 in y2s
                if norm([x2, y2] - x) <= dist
                    push!(cps, [x2, y2])
                    push!(dists, norm([x2, y2] - x))
                end
            end
        end

        return cps[sortperm(dists)]
    end

    if length(x) == 3
        x2s = Int.(floor(x[1] - dist) : 1 : ceil(x[1] + dist))
        y2s = Int.(floor(x[2] - dist) : 1 : ceil(x[2] + dist))
        z2s = Int.(floor(x[3] - dist) : 1 : ceil(x[3] + dist))

        for x2 in x2s
            for y2 in y2s
                for z2 in z2s
                    if norm([x2, y2, z2] - x) <= dist
                        push!(cps, [x2, y2, z2])
                        push!(dists, norm([x2, y2, z2] - x))
                    end
                end
            end
        end

        return cps[sortperm(dists)]
    end
end

## Note that we assume that there is no degeneracy because of randomness
r = 5 # radius
num_samples = 1e3 # number of samples
num_points = 20
for n = [2, 3] # dimension
    for _ in 1 : num_samples
        
        # In order to find certain number of closest points, 
        # the function `all_closest_points_within_distance` brute force search a large enough region, 
        # followed by sorting the results
        #
        # Because of randomness, no two points will have the same distance from the input vector
        # hence the ordering of the lattice points is unique
        x = r * rand(n) .- r/2
        cps = all_closest_points_within_distance(x, r)[1:num_points]
        cps2 = closest_points_Zn(x, num_points)
        @test cps == cps2
    end
end

# For input vectors that are integer valued or half-integer valued,
# there are cases where two lattice points share the same distance
# Here we test that closest_points_Zn return no repeated result,
# and the returned lattice points have non-decreasing distance.
test_cases = [[0, 0], [0, 1], [1/2, 0], [1/2, 1/2, 1/2], [0, 1/2, 1]]

for x in test_cases
    cps = closest_points_Zn(x, num_points)
    @test length(cps) == length(Set(cps)) # no repetition

    dists = [norm(cp - x) for cp in cps]
    @test (dists[2:end] .- dists[1:end-1] .>= 0) == ones(length(dists)-1) # dist is non-decreasing
end

## Tests for closest_points_Dn_dual for integer and half-integer inputs

M = Dn_dual(2)
cps = closest_points_Dn_dual([0.5, 0.5], 9)
@test Set(cps) == Set([[0.5, 0.5], [0, 0], [0, 1], [1, 0], [1, 1], [0.5, 1.5], [1.5, 0.5], [0.5, -0.5], [-0.5, 0.5]])
cps = closest_points_Dn_dual([0., 0.], 9)
@test Set(cps) == Set([[0, 0], [-0.5, -0.5], [0.5, -0.5], [-0.5, 0.5], [0.5, 0.5], [1, 0], [0, 1], [-1, 0], [0, -1]])
cps = closest_points_Dn_dual([0.5, -0.5], 9)
@test Set(cps) == Set([[0.5, -0.5], [0, -1], [1, -1], [0, 0], [1, 0], [1.5, -0.5], [0.5, 0.5], [-0.5, -0.5], [0.5, -1.5]])

M = Dn_dual(3)
cps = closest_points_Dn_dual([0.5, 0.5, 0.5], 9)
@test Set(cps) == Set([[0.5, 0.5, 0.5], [0, 0, 0], [1, 0, 0], [0, 0, 1], [1, 0, 1], [0, 1, 0], [1, 1, 1], [0, 1, 1], [1, 1, 0]])
cps = closest_points_Dn_dual([0., 0., 0.], 9)
@test Set(cps) == Set([[0, 0, 0], [-0.5, -0.5, -0.5], [0.5, -0.5, -0.5], [-0.5, -0.5, 0.5], [0.5, -0.5, 0.5], [-0.5, 0.5, -0.5], [0.5, 0.5, 0.5], [-0.5, 0.5, 0.5], [0.5, 0.5, -0.5]])
cps = closest_points_Dn_dual([0.5, 0.5, -0.5], 9)
@test Set(cps) == Set([[0.5, 0.5, -0.5], [0, 0, -1], [1, 0, -1], [0, 0, 0], [1, 0, 0], [0, 1, -1], [1, 1, 0], [0, 1, 0], [1, 1, -1]])



test_cases = [[0, 0], [0, 1], [1/2, 0], [1/2, 1/2, 1/2], [0, 1/2, 1]]
num_points = 20
for x in test_cases
    cps = closest_points_Dn_dual(x, num_points)
    @test length(cps) == length(Set(cps)) # no repetition

    dists = [norm(cp - x) for cp in cps]
    @test (dists[2:end] .- dists[1:end-1] .>= 0) == ones(length(dists)-1) # dist is non-decreasing
end

## Tests for closest_points_Dn_dual for random inputs

r = 5 # radius
num_samples = 1e3 # number of samples
num_points = 20
for n = [2, 3] # dimension
    for _ in 1 : num_samples
        x = r * rand(2) .- r/2
        cps = all_closest_points_within_distance(x, r)[1:num_points]
        r1 = 0.5 * ones(length(x))
        cps2 = all_closest_points_within_distance(x-r1, r)[1:num_points]
        cps2 = [cp + r1 for cp in cps2]
        
        cps3 = vcat(cps, cps2)
        dists3 = [norm(cp-x) for cp in cps3]
        
        cps3 = cps3[sortperm(dists3)][1:num_points]

        cps4 = closest_points_Dn_dual(x, num_points)
        @test cps3 == cps4
    end
end



## Tests for closest_points_Dn_dual for integer and half-integer inputs

test_cases = [[0, 0], [0, 1], [1/2, 0], [1/2, 1/2, 1/2], [0, 1/2, 1]]
num_points = 20
for x in test_cases
    cps = closest_points_Dn(x, num_points)
    @test length(cps) == length(Set(cps)) # no repetition

    dists = [norm(cp - x) for cp in cps]
    @test (dists[2:end] .- dists[1:end-1] .>= 0) == ones(length(dists)-1) # dist is non-decreasing
    for cp in cps
        @test mod(sum(cp), 2) == 0
    end
end

## Tests for closest_points_Dn_dual for random inputs
r = 5 # radius
num_samples = 1e3 # number of samples
num_points = 20
for n = [2, 3] # dimension
    for _ in 1 : num_samples
        
        x = r * rand(2) .- r/2
        cps = all_closest_points_within_distance(x, r)

        cps2 = []
        for cp in cps
            if mod(sum(cp), 2) == 0
                push!(cps2, cp)
            end
        end
        cps2 = cps2[1:num_points]


        cps3 = closest_points_Dn(x, num_points)
        @test cps2 == cps3
    end
end