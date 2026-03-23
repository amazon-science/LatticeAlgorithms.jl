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
using DataStructures
using LinearAlgebra
using SparseArrays

# Test example from https://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-greedy-algo-7/
g = [[0, 4, 0, 0, 0, 0, 0, 8, 0],
    [4, 0, 8, 0, 0, 0, 0, 11, 0],
    [0, 8, 0, 7, 0, 4, 0, 0, 2],
    [0, 0, 7, 0, 9, 14, 0, 0, 0],
    [0, 0, 0, 9, 0, 10, 0, 0, 0],
    [0, 0, 4, 14, 10, 0, 2, 0, 0],
    [0, 0, 0, 0, 0, 2, 0, 1, 6],
    [8, 11, 0, 0, 0, 0, 1, 0, 7],
    [0, 0, 2, 0, 0, 0, 6, 7, 0]
    ]

g = Matrix(hcat(g...))
result = []
for src2 in 2:size(g, 1)
    push!(result, shortest_path(g, 1, src2))
end
@test [4, 12, 19, 21, 11, 9, 8, 14] == [res[1] for res in result]

for i in 1 : length(result)
    dist = result[i][1]
    path = result[i][2]
    dist_path = sum([g[path[j], path[j+1]] for j in 1 : length(path)-1])
    @test dist == dist_path

    # path2 = result2[i][2]
    # dist_path2 = sum([g2[e][2] for e in path2])
    # @test dist == dist_path2
end

# Test examples from https://github.com/JuliaGraphs/Graphs.jl/blob/master/test/shortestpaths/dijkstra.jl
w = [
    0.0 3.0 0.0 1.0
    3.0 0.0 2.0 0.0
    0.0 2.0 0.0 3.0
    1.0 0.0 3.0 0.0
]

result = [shortest_path(w, i, 2) for i in [1,3,4]]
@test [res[2] for res in result] == [[1, 2], [3, 2], [4, 1, 2]]
result2 = shortest_paths(w, 2, 4, 2)
@test result2[1] == [4.0, 5.0]
@test result2[2] == [[2, 1, 4], [2, 3, 4]]
# # Test the new data structure of g
# w = get_graph_from_matrix(w)
# result_new = [shortest_path(w, i, 2) for i in [1,3,4]]


w2 = [
    0.0 3.0 0.0 1.0
    3.0 8.0 2.0 0.0
    0.0 2.0 0.0 3.0
    1.0 0.0 3.0 0.0
]
result3 = [shortest_path(w2, i, 2) for i in [1,3,4]]
@test result3 == result
result4 = shortest_paths(w2, 2, 4, 3)
@test result4[1] == [4.0, 5.0]
@test result4[2] == [[2, 1, 4], [2, 3, 4]] # There are only two shortest paths


w3 = [
    0.0 3.0 0.0 1.0
    3.0 8.0 2.0 1.0
    0.0 2.0 0.0 3.0
    1.0 1.0 3.0 0.0
]
result5 = shortest_paths(w3, 2, 4, 3)
@test result5[1] == [1.0, 4.0, 5.0]
@test result5[2] == [[2, 4], [2, 1, 4], [2, 3, 4]]

# Testing with the example in https://en.wikipedia.org/wiki/Yen%27s_algorithm
w4 = [
    0.0 3.0 2.0 0.0 0.0 0.0
    0.0 0.0 0.0 4.0 0.0 0.0
    0.0 1.0 0.0 2.0 3.0 0.0
    0.0 0.0 0.0 0.0 2.0 1.0
    0.0 0.0 0.0 0.0 0.0 2.0
    0.0 0.0 0.0 0.0 0.0 0.0
]
result6 = shortest_paths(w4, 1, 6, 3, true)
@test result6[1] == [5.0, 7.0, 8.0]
@test result6[2][1] == [1, 3, 4, 6]
@test result6[2][2] == [1, 3, 5, 6]
@test result6[2][3] ∈ [[1, 3, 4, 5, 6], [1, 3, 2, 4, 6], [1, 2, 4, 6]]

w5 = [
    0.0 3.0 2.0 3.0 0.0 0.0
    0.0 0.0 0.0 4.0 0.0 0.0
    0.0 1.0 0.0 2.0 3.0 0.0
    3.0 0.0 0.0 0.0 2.0 1.0
    0.0 0.0 0.0 0.0 0.0 2.0
    0.0 0.0 0.0 0.0 0.0 0.0
]
result7 = shortest_paths(w5, 1, 6, 100, true)
@test result7[1] == [4.0, 5.0, 7.0, 7.0, 8.0, 8.0, 8.0, 11.0, 11.0]

m = float([0 2 2 0 0 1; 2 0 1 0 0 0; 2 1 0 4 0 0; 0 0 4 0 1 0; 0 0 0 1 0 1; 1 0 0 0 1 0])
result = [shortest_path(m, 3, i) for i in 1:size(m,1)]
dists = [res[1] for res in result]
@test dists[[1,2,3,6]] == [2,1,0,3]
@test dists[4] > 3
@test dists[5] > 3


w0 = [
    0.0 3.0 0.0 0.0
    3.0 8.0 0.0 0.0
    0.0 0.0 0.0 3.0
    0.0 0.0 3.0 0.0
]
@test shortest_path(w0, 1, 4) == (Inf, [])
# # Test the new data structure of g
# w0 = get_graph_from_matrix(w0)
# @test shortest_path(w0, 1, 4) == (Inf, [])

## More tests for functions for finding the shortest paths

"""
    brute_force_search_shortest_paths(g::Matrix, source::Int, target::Int, directed_graph::Bool=false)
    
Brute force search all the paths from source to target, then sort them by weights
"""
function brute_force_search_shortest_paths(g::Matrix, source::Int, target::Int, directed_graph::Bool=false)
    if source == target
        dists = [0]
        paths = [[source]]
        return dists, paths
    end

    paths_to_explore = [[source]]

    dists = []
    paths = []
    while length(paths_to_explore) > 0
        path = pop!(paths_to_explore)
        last_element = path[end]
        if directed_graph
            outneighbors2 = findall(g[last_element,last_element+1:end] .> 0)
        else
            outneighbors2 = findall(g[last_element,:] .> 0)
        end
        for v in outneighbors2
            path2 = vcat(path, [v])
            if v == target
                dist = sum([g[path2[i], path2[i+1]] for i in 1 : length(path2)-1])
                push!(dists, dist)
                push!(paths, path2)
            elseif length(path2) == length(Set(path2))
                push!(paths_to_explore, path2)
            end
        end
    end
    sortind = sortperm(dists)
    dists = dists[sortind]
    paths = paths[sortind]
    return dists, paths
end

test_cases = [[2, 5], [1, 3], [4, 8], [7, 6]]

num_samples = 100
for _ in num_samples
    g = rand(8, 8) * 4
    g = g + transpose(g)
    g = g - diagm(diag(g))

    for case in test_cases
        source, target = case[1], case[2]
        res1 = brute_force_search_shortest_paths(g, source, target)
        res2 = shortest_paths(g, source, target, length(res1[2]))
        # @test res1[1] == res2[1]
        @test all((res1[1] - res2[1]) .< 1e-10)
        @test all([r ∈ res2[2] for r in res1[2]])
    end
end


## Tests for functions for finding minimum weight cycles

# Test example from https://www.geeksforgeeks.org/find-minimum-weight-cycle-undirected-graph/

g = [[0, 4, 0, 0, 0, 0, 0, 8, 0],
    [4, 0, 8, 0, 0, 0, 0, 11, 0],
    [0, 8, 0, 7, 0, 4, 0, 0, 2],
    [0, 0, 7, 0, 9, 14, 0, 0, 0],
    [0, 0, 0, 9, 0, 10, 0, 0, 0],
    [0, 0, 4, 14, 10, 0, 2, 0, 0],
    [0, 0, 0, 0, 0, 2, 0, 1, 6],
    [8, 11, 0, 0, 0, 0, 1, 0, 7],
    [0, 0, 2, 0, 0, 0, 6, 7, 0]
    ]
g = Matrix(hcat(g...))
g = Float64.(g)

c1 = minimum_weight_cycle(g)
@test c1[1] == 14

num_samples = 100
for _ in 1:num_samples
    g = rand(8, 8) * 4
    g = g + transpose(g)
    g = g - diagm(diag(g))

    res_list = []
    for i in 1 : size(g, 1)
        for j in i+1:size(g, 1)
            res = brute_force_search_shortest_paths(g, i, j)
            dists, paths = res[1], res[2]
            dist = dists[1] + dists[2]
            cycle = vcat(paths[1], paths[2])
            push!(res_list, ((i, j), dist, paths[1], paths[2]))
        end
    end
    dists = [r[2] for r in res_list]
    paths = [(r[3], r[4]) for r in res_list]

    c1 = minimum_weight_cycle(g)
    cycle1 = [[c1[2][j], c1[2][j+1]] for j in 1 : length(c1[2])-1]

    ind = argmin(dists)

    @test c1[1] ≈ dists[ind]

    cycle2 = []
    for path in paths[ind]
        for j in 1 : length(path)-1
            push!(cycle2, [path[j], path[j+1]])
        end
    end

    cycle1 = sort(sort.(cycle1))
    cycle2 = sort(sort.(cycle2))
    if !(cycle1 == cycle2)
        println(g)
    end
    @test cycle1 == cycle2
end
