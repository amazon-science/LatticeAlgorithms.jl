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

M = [1 -1; 1 1]

@test closest_point([0, 0.5], M) == [0.0, 0.0]
@test norm(shortest_vector(M)) ≈ √2


# massive tests

function test_closest_point(N::Integer, option::String, num_samples=1e2)
    elapsed_time = @elapsed for _ in 1 : num_samples

        M = rand(N, N) * N .- N/2

        if option == "LLL"
            basis = lll(M)
        elseif option == "KZ"
            basis = kz(M)
        elseif option == "LQ"
            basis = lq_reduce(M)
        else
            basis = M
        end
        
        xs = [rand(N) * N .- N/2 for _ in 1 : num_samples]
        us = [closest_point(x, basis) for x in xs]
        ds = [norm(u-x) for (u, x) in zip(us, xs)]
        ds2 = [
            findmin(
                vcat(
                    [norm(M[i, :] - x) for i in 1:N], 
                    [norm(M[i, :] + x) for i in 1:N]
                )
            )[1]
            for x in xs
        ]
        # to avoid floating errors
        ds = round.(ds, digits = 10)
        ds2 = round.(ds2, digits = 10) 

        @test sum(ds .<= ds2) == num_samples
    end
    # println("test_closest_point_$option for N = $N, elapsed_time=$elapsed_time")
end

for N in 2 : 5
    test_closest_point(N, "LLL")
    test_closest_point(N, "KZ")
    test_closest_point(N, "LQ")
    test_closest_point(N, "NONE")
end


function test_shortest_vector(N::Integer, num_samples=1e2)
    elapsed_time = @elapsed for _ in 1 : num_samples

        M = rand(N, N) * N .- N/2
        sv = shortest_vector(M)
        d2 = findmin([norm(M[i, :]) for i in 1:N])[1]

        # to avoid floating errors
        d1 = round(norm(sv), digits=10)
        d2 = round(d2, digits = 10)

        @test d1 <= d2

    end
    # println("test_shortest_vector$option for N = $N, elapsed_time=$elapsed_time")
end

for N in 2 : 5
    test_shortest_vector(N)
end

