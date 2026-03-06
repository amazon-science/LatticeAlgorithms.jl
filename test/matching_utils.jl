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


## Tests for functions for finding shortest paths

# function get_graph_from_matrix(w0)
#     g = Vector{Tuple{Vector{Int64}, Float64}}()
#     for i in 1 : size(w0, 1)
#         for j in i+1 : size(w0, 2)
#             w0[i, j] > 0 && push!(g, ([i, j], w0[i, j]))
#         end
#     end
#     return g
# end

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
# # Test the different data structure of g
# g2 = get_graph_from_matrix(g)
# result2 = []
# for src2 in 2:size(g, 1)
#     push!(result2, shortest_path(g2, 1, src2))
# end
# @test [4, 12, 19, 21, 11, 9, 8, 14] == [res[1] for res in result2]


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
# # Test the new data structure of g
# w2 = get_graph_from_matrix(w2)
# result3_new = [shortest_path(w2, i, 2) for i in [1,3,4]]
# @test result3_new == result_new


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
# # Test the new data structure of g
# m2 = get_graph_from_matrix(m)
# result = [shortest_path(m2, 3, i) for i in 1:size(m,1)]
# dists = [res[1] for res in result]
# @test dists[[1,2,3,6]] == [2,1,0,3]
# @test dists[4] > 3
# @test dists[5] > 3

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
# # Test the new data structure of g
# g2 = get_graph_from_matrix(g)
# c1 = minimum_weight_cycle(g2)
# @test c1[1] == 14


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

    # # Test the new data structure of g
    # g2 = get_graph_from_matrix(g)
    # c3 = minimum_weight_cycle(g2)[2]
    # cycle3 = [g2[e][1] for e in c3]
    # cycle3 = sort(sort.(cycle3))
    # @test cycle3 == cycle2
end

# Test mwms, making sure that the weights of the matchings are ordered
σ = 0.6
K = 20
num_samples = 1000
for d in [5, 7]
    stabilizers = surface_code_Z_stabilizers(d)
    for ind_sample in 1 : num_samples
        ηs = σ * randn(d^2)
        g, hvs, _, _, _, G, vertices_edge_mapping = decoding_graph(ηs, stabilizers)
        Matchings, _ = mwms(g, hvs, G, vertices_edge_mapping, K)
        weights = [length(M1)==0 ? 0 : sum([g[ee[1], ee[2]] for ee in M1]) for M1 in Matchings]
        @assert issorted(weights) == true
    end
end

# ## Testing mwm with the new data structure for g
# num_samples = 1e3
# for d in [3, 5, 7]
#     stabilizers = surface_code_X_stabilizers(d)
#     for _ in 1 : num_samples
        
#         x = 0.6 * randn(d^2)
#         η = 1

#         # get the weights for each qubit, which are the edges of the graph
#         closest_integers = closest_integer.(x)
#         second_closest_integers = second_closest_integer.(x)

#         edge_weight_list = ((second_closest_integers.-x).^2 - (closest_integers.-x).^2) .* η.^2

#         # get edge_to_mode_dict
#         num_vertices = length(stabilizers)+1
#         keys_stabilizers = collect(keys(stabilizers))
#         values_stabilizers = collect(values(stabilizers))
#         edge_to_mode_dict = Dict()
#         for qubit in 1 : length(x)
#             qubit_in_stab = findall(qubit .∈ values_stabilizers) # determine which stabilizers the qubit is in
#             if length(qubit_in_stab) > 2
#                 error("Cannot decode code where a single fault can lead more than 2 errors.")
#             elseif length(qubit_in_stab) == 2
#                 vertex_1, vertex_2 = keys_stabilizers[qubit_in_stab[1]], keys_stabilizers[qubit_in_stab[2]]
#             elseif length(qubit_in_stab) == 1
#                 vertex_1, vertex_2 = keys_stabilizers[qubit_in_stab[1]], num_vertices
#             elseif length(qubit_in_stab) == 0 # ignore if the qubit is not in any stabilizers
#                 continue
#             end

#             if Set([vertex_1, vertex_2]) in keys(edge_to_mode_dict)
#                 if edge_weight_list[qubit] < edge_to_mode_dict[Set([vertex_1, vertex_2])][2]
#                     edge_to_mode_dict[Set([vertex_1, vertex_2])] = (qubit, edge_weight_list[qubit])
#                 end
#             else
#                 merge!(edge_to_mode_dict, Dict(Set([vertex_1, vertex_2]) => (qubit, edge_weight_list[qubit])))
#             end
#         end

#         # highlight the unhappy stabilizers/vertices
#         highlighted_vertices = zeros(Int, num_vertices)
#         for (index, stabilizer) in stabilizers
#             if mod(sum(closest_integers[stabilizer]), 2) == 1
#                 highlighted_vertices[index] = 1
#             end
#         end

#         highlighted_vertices[num_vertices] = mod(sum(highlighted_vertices), 2)

#         # Get the weights of the graph
#         g = Matrix(0.0I, num_vertices, num_vertices)
#         for (key, (qubit, weight)) in edge_to_mode_dict
#             vertex_1, vertex_2 = collect(key)[1], collect(key)[2]
#             g[vertex_1, vertex_2] = weight
#         end
#         g = g + transpose(g)

#         g2 = get_graph_from_matrix(g)

#         M1 = mwpm(g, highlighted_vertices)    
#         M1 = sort(sort.(collect.(M1)))

#         M2 = mwm(g2, highlighted_vertices)
#         M2 = [g2[idx][1] for idx in M2]
#         length(M1) > 0 && (@test M1 == M2)
#     end
# end


# ## Tests using mwms to find multiple minimum weight cycles

# function combinations(arr, k)
#     k == 1 && return [[i] for i in arr]
#     k == length(arr) && return [arr]
    
#     result = []
#     for i in 1 : length(arr)-1
#         res2 = combinations(arr[i+1:end], k-1)
#         res2 = [vcat([arr[i]], res) for res in res2]
#         result = vcat(result, res2)
#     end
#     return result
# end

# @test combinations([1,2], 2) == [[1, 2]]

# @test combinations([1,2,3], 1) == [[1], [2], [3]]

# @test combinations([1,2,3], 2) == [[1, 2], [1, 3], [2, 3]]

# @test combinations([1,2,3,4,5], 3) == [
#     [1, 2, 3], 
#     [1, 2, 4], 
#     [1, 2, 5], 
#     [1, 3, 4], 
#     [1, 3, 5], 
#     [1, 4, 5], 
#     [2, 3, 4], 
#     [2, 3, 5], 
#     [2, 4, 5], 
#     [3, 4, 5]
# ]

# function brute_force_search_minimum_weight_cycles(g::Matrix, K::Int)
#     cycles_list = []
#     edges = findall(triu(g, 1).!=0) # exclude the main diagonal
#     for edge in edges
#         g2 = deepcopy(g)

#         weight_edge = g2[edge[1], edge[2]]
#         g2[edge[1], edge[2]] = 0
#         g2[edge[2], edge[1]] = 0

#         _, paths = shortest_paths(g2, edge[1], edge[2], K)
#         paths = [vcat(path, [path[1]]) for path in paths]

#         cycles = [sort(sort.([[c[j], c[j+1]] for j in 1 : length(c)-1])) for c in paths]

#         cycles_list = vcat(cycles_list, cycles)
#     end
#     cycles = unique(cycles_list)
#     weights = [sum([g[e[1], e[2]] for e in cycle]) for cycle in cycles]
#     ind = sortperm(weights)
#     weights = weights[ind]
#     cycles = cycles[ind]

#     if length(weights) > K
#         weights, cycles = weights[1:K], cycles[1:K]
#     end

#     # Up to now, c1, c2 ∈ cycles are not subset of each other if c1 ≠ c2.
#     # Now we need to consider multiple cycles, such as c3 = c1 ∪ c2
#     # if w(c3) is smaller than the weight of the last cycle in cycles,
#     # then c3 should be placed somewhere in the cycles

#     # Find the largest index `ind` such that combining cycles[1] and cycles[ind]
#     # yield a cycle with weight less than weights[end]
#     # We note that w(cycles[i] ⊕ cycles[j]) <= w(cycles[i]) + w(cycles[j])
#     # because edges happens even number of times cancel each other

#     ind = 0 
#     for (ind2, cycle) in enumerate(cycles[2:end])
#         # println(ind2)
#         matchings = [cycles[1], cycle]
#         matching = vcat(matchings...)
#         matching = sort(sort.(matching))
#         matching_counter = counter(matching)
#         matching = Vector{Vector{Int64}}()
#         for (key, val) in matching_counter
#             if mod(val, 2)≠0
#                 push!(matching, key)
#             end
#         end
#         matching = sort(sort.(matching))
#         w_matching = sum([g[e[1], e[2]] for e in matching])
#         if matching ∉ cycles && w_matching < weights[end]
#             ind = ind2+1
#         end
#     end

#     # ind = findall(cumsum(weights) .> weights[end])[1]-1
    
#     ind == 1 && return weights, cycles

#     new_cycles =[ ]
#     new_weights =[ ]
#     for num in 2 : ind # number of matchings
#         # println(num)
#         combs = combinations(collect(1:ind), num) # combinations with num of matchings that has weight less than weights2[end]
#         for comb in combs
#             # println(comb)
#             matchings = cycles[comb]
#             matching = vcat(matchings...)
#             matching = sort(sort.(matching))
#             matching_counter = counter(matching)
#             matching = Vector{Vector{Int64}}()
#             for (key, val) in matching_counter
#                 if mod(val, 2)≠0
#                     push!(matching, key)
#                 end
#             end
#             if length(matching) > 0 
#                 matching = sort(sort.(matching))
#                 w_matching = sum([g[e[1], e[2]] for e in matching])
#                 if matching ∉ cycles && matching ∉ new_cycles && w_matching < weights[end]
#                     # println(matching)
#                     push!(new_cycles, matching)
#                     push!(new_weights, w_matching)
#                 end
#             end
#         end
#     end
#     weights = vcat(weights, new_weights)
#     cycles = vcat(cycles, new_cycles)
#     ind2 = sortperm(weights)
#     weights = weights[ind2]
#     cycles = cycles[ind2]
#     if length(weights) > K
#         weights, cycles = weights[1:K], cycles[1:K]
#     end

#     return weights, cycles
# end

# num_samples = 100
# Nrange = [8, 10, 12]
# num_MWCs = 20
# for N in Nrange
#     println("N, num_samples, num_MWCs = $N, $(num_samples), $(num_MWCs)")
#     # flush(stdout)
#     for iii in 1:num_samples 
#         # println(iii)
#         # flush(stdout)
#         g = rand(N, N) * N/2
#         g = g + transpose(g)
#         g = g - diagm(diag(g))
        
#         syndrome = zeros(Int, size(g, 1))

#         weights, cycles = brute_force_search_minimum_weight_cycles(g, num_MWCs)

#         cycles2 = mwms(g, syndrome, num_MWCs)

#         weights2 = [sum([g[e[1], e[2]] for e in cycle]) for cycle in cycles2]
#         if !(all(weights .≈ weights2))
#             println(g)
#         end
#         @test all(weights .≈ weights2)
#         @test cycles == cycles2 # assume no degeneracy
        
#         diff_weights = [weights[i+1] - weights[i] for i in 1 : length(weights)-1]
#         @test all(diff_weights .>= 0)
#         @test length(Set(cycles2)) == length(cycles2) # No repeated matching

        
#         # Testing secondmwm with the old and new data structures of g
#         c2 = secondmwm(g, syndrome)[2]
#         # g3 = get_graph_from_matrix(g)
#         # c3 = secondmwm(g3, syndrome)[2]
#         # c3 = [g3[e][1] for e in c3]
#         # println(c3)
#         # c3 = sort.(sort(c3))
#         @test cycles2[2] == c2
#         # @test cycles2[2] == c3

#         # # # Testing the new data structure for g
        
#         # cycles3 = mwms(g3, syndrome, num_MWCs)
#         # cycles3 = [[g3[e][1] for e in c] for c in cycles3]
#         # @test cycles2 == cycles3 # assume no degeneracy
#     end
# end


# ## secondmwm
# g = [[0, 4, 0, 0, 0, 0, 0, 8, 0],
#     [4, 0, 8, 0, 0, 0, 0, 11, 0],
#     [0, 8, 0, 7, 0, 4, 0, 0, 2],
#     [0, 0, 7, 0, 9, 14, 0, 0, 0],
#     [0, 0, 0, 9, 0, 10, 0, 0, 0],
#     [0, 0, 4, 14, 10, 0, 2, 0, 0],
#     [0, 0, 0, 0, 0, 2, 0, 1, 6],
#     [8, 11, 0, 0, 0, 0, 1, 0, 7],
#     [0, 0, 2, 0, 0, 0, 6, 7, 0]
#     ]

# g = Matrix(hcat(g...))

# syndrome = [1, 0, 0, 0, 0, 0, 0, 0, 1]
# w_M2, M2 = secondmwm(g, syndrome)
# @test M2 == sort(sort.([[8, 9], [1, 8]]))
# @test w_M2 == 15 # weight of M1 == 14, hence weight(M2)>=15
# # Test the new data structure for g
# g2 = get_graph_from_matrix(g)
# w_M2_2, M2_2 = secondmwm(g2, syndrome)
# M2_2 = [g2[e][1] for e in M2_2]
# @test sort(sort.(M2_2)) == sort(sort.([[8, 9], [1, 8]]))
# @test w_M2_2 == 15

# w = [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0]
# syndrome = [1,0,0,1]
# w_M2, M2 = secondmwm(w, syndrome)
# @test M2 == sort(sort.([[4, 3], [3, 2], [2, 1]]))
# @test w_M2 == 3
# # Test the new data structure for g
# g2 = get_graph_from_matrix(w)
# w_M2_2, M2_2 = secondmwm(g2, syndrome)
# M2_2 = [g2[e][1] for e in M2_2]
# @test sort(sort.(M2_2)) == sort(sort.([[4, 3], [3, 2], [2, 1]]))
# @test w_M2_2 == 3


# # Test example from https://www.geeksforgeeks.org/difference-between-the-shortest-and-second-shortest-path-in-an-unweighted-bidirectional-graph/
# w = [
#     [0, 1, 0, 1],
#     [0, 0, 1, 0],
#     [0, 0, 0, 1],
#     [0, 0, 0, 0]]

# w = Matrix(hcat(w...))
# w = w + transpose(w)
# syndrome = [1,0,0,1]
# w_M2, M2 = secondmwm(w, syndrome)
# @test w_M2 == 3
# @test M2 == sort(sort.([[4, 3], [3, 2], [2, 1]]))
# # Test the new data structure for g
# g2 = get_graph_from_matrix(w)
# w_M2_2, M2_2 = secondmwm(g2, syndrome)
# M2_2 = [g2[e][1] for e in M2_2]
# @test sort(sort.(M2_2)) == sort(sort.([[4, 3], [3, 2], [2, 1]]))
# @test w_M2_2 == 3


# # Test example from https://www.geeksforgeeks.org/path-from-a-given-source-to-a-given-destination-having-kth-largest-weight-in-a-graph/
# w = [
#     [0, 10, 0 , 40, 0, 0, 0],
#     [0,  0, 10,  0, 0, 0, 0],
#     [0,  0,  0, 10, 0, 0, 0],
#     [0,  0,  0,  0, 2, 0, 0],
#     [0,  0,  0,  0, 0, 3, 8],
#     [0,  0,  0,  0, 0, 0, 3],
#     [0,  0,  0,  0, 0, 0, 0]]

# w = Matrix(hcat(w...))
# w = w + transpose(w)
# syndrome = [1,0,0,0,0,0,1]
# w_M2, M2 = secondmwm(w, syndrome)
# @test w_M2 == 40
# @test M2 == sort(sort.([[5, 7], [2, 3], [1, 2], [4, 5], [3, 4]]))
# # Test the new data structure for g
# g2 = get_graph_from_matrix(w)
# w_M2_2, M2_2 = secondmwm(g2, syndrome)
# M2_2 = [g2[e][1] for e in M2_2]
# @test sort(sort.(M2_2)) == sort(sort.([[5, 7], [2, 3], [1, 2], [4, 5], [3, 4]]))
# @test w_M2_2 == 40


# ## Test secondmwm for arbitrary graphs with only two syndromes

# test_cases = [[2, 5], [1, 3], [4, 8], [7, 6]]
# num_samples = 1000
# N = 8 
# for _ in 1:num_samples
#     g = rand(N, N) * N/2
#     g = g + transpose(g)
#     g = g - diagm(diag(g))

#     for case in test_cases
#         # println(case)
#         source, target = case[1], case[2]
#         res1 = shortest_paths(g, source, target, 2)

#         wc, c = minimum_weight_cycle(g)

#         if res1[1][2] > wc + res1[1][1] 
#             # the second mwm is the shortest path with the mwc
#             M2_v2 = res1[2][1]
#             M2_v2 = [[M2_v2[j], M2_v2[j+1]] for j in 1 : length(M2_v2)-1]
            
#             M2_v2 = [M2_v2..., [[c[j], c[j+1]] for j in 1 : length(c)-1]...]
            
#         else
#             # the second mwm is the 2nd shortest path
#             M2_v2 = res1[2][2]
#             M2_v2 = [[M2_v2[j], M2_v2[j+1]] for j in 1 : length(M2_v2)-1]
#         end
#         M2_v2 = sort(sort.(M2_v2))
#         w_M2_v2 = sum([g[e[1], e[2]] for e in M2_v2])

#         syndrome = zeros(Int, size(g, 1))
#         syndrome[source] = syndrome[target] = 1
#         w_M2, M2 = secondmwm(g, syndrome)

#         if !(w_M2 ≈ w_M2_v2)
#             println(g)
#             println(syndrome)
#         end        
#         @test w_M2 ≈ w_M2_v2
#         @test M2 == M2_v2

#         # Test the new data structure for g
#         g2 = get_graph_from_matrix(g)
#         w_M2_2, M2_2 = secondmwm(g2, syndrome)
#         M2_2 = sort(sort.([g2[e][1] for e in M2_2]))
#         @test w_M2 ≈ w_M2_2
#         @test M2 == M2_2
#     end
# end

# ## Test secondmwm and mwms (K=2) for arbitrary graphs with more syndromes

# # First define a function that partition the syndromes into multiple pairs
# function partition_array(arr)
#     length(arr) == 2 && return [[arr]]
    
#     subcases = []
#     for i in 2 : length(arr)
#         subcase1 = [arr[1], arr[i]]
#         # println("subcase1 = $subcase1")
#         subcase2 = setdiff(arr, subcase1)
        
#         subcase3 = partition_array(subcase2)
        
#         # println("subcase3 = $subcase3")
#         for subcase in subcase3
#             # println("subcase = $subcase")
            
#             # push!(subcases, push!(subcase, subcase1))
#             push!(subcases, [subcase..., subcase1])
#         end        
#     end
#     # println("subcases = $subcases")
#     # println()
#     return subcases
# end
# @test partition_array([1,2]) == [[[1,2]]]
# @test partition_array([1,2,3,4]) == [
#     [[3, 4], [1, 2]],
#     [[2, 4], [1, 3]],
#     [[2, 3], [1, 4]],
# ]

# @test partition_array([1,2,3,4,5,6]) == [
#     [[5, 6], [3, 4], [1, 2]],
#     [[4, 6], [3, 5], [1, 2]],
#     [[4, 5], [3, 6], [1, 2]],
#     [[5, 6], [2, 4], [1, 3]],
#     [[4, 6], [2, 5], [1, 3]],
#     [[4, 5], [2, 6], [1, 3]],
#     [[5, 6], [2, 3], [1, 4]],
#     [[3, 6], [2, 5], [1, 4]],
#     [[3, 5], [2, 6], [1, 4]],
#     [[4, 6], [2, 3], [1, 5]],
#     [[3, 6], [2, 4], [1, 5]],
#     [[3, 4], [2, 6], [1, 5]],
#     [[4, 5], [2, 3], [1, 6]],
#     [[3, 5], [2, 4], [1, 6]],
#     [[3, 4], [2, 5], [1, 6]],
# ]

# num_samples = 1000
# Nrange = [8, 10, 12]
# # Nrange = [8, 10, 12, 14, 16, 18] # Tested but take a while. all passed

# for N in Nrange
#     println("N=$N")
#     for _ in 1:num_samples 
#         g = rand(N, N) * N/2
#         g = g + transpose(g)
#         g = g - diagm(diag(g))

#         syndrome = round.(Int, rand(size(g, 1)-1))
#         if mod(sum(syndrome), 2) == 1
#             push!(syndrome, 1)
#         else
#             push!(syndrome, 0)
#         end
        
#         syndrome == zeros(Int, length(syndrome)) && continue

#         w_M2_true, M2_true = secondmwm(g, syndrome)
#         case = findall(syndrome .== 1)

#         # Test the new data structure for g
#         g2 = get_graph_from_matrix(g)
#         w_M2_2, M2_2 = secondmwm(g2, syndrome)
#         M2_2 = sort(sort.([g2[e][1] for e in M2_2]))
#         @test w_M2_true ≈ w_M2_2
#         @test M2_true == M2_2

#         subcases = partition_array(case)

#         paths_list2 = []
#         for subcase in subcases
#             dists_list = []
#             paths_list = []
#             for subsubcase in subcase
#                 source, target = subsubcase[1], subsubcase[2]
#                 dists, paths = shortest_paths(g, source, target, 2)
#                 paths = [[[path[j], path[j+1]] for j in 1 : length(path)-1] for path in paths]
#                 push!(dists_list, dists)
#                 push!(paths_list, paths)
#                 # println(dist)
#                 # println(path)
#                 # println()
#             end
        
#             M2_1 = vcat([paths[1] for paths in paths_list]...)
#             M2_1 = sort(sort.(M2_1))
#             M2_1_counter = counter(M2_1)
#             M2_1 = []
#             for (key, val) in M2_1_counter
#                 if mod(val, 2)≠0
#                     push!(M2_1, key)
#                 end
#             end
            
#             diff_list = [dists[2] - dists[1] for dists in dists_list]
#             diff_min, diff_min_ind = findmin(diff_list)
#             M2_2 = vcat([ind == diff_min_ind ? paths[2] : paths[1] for (ind, paths) in enumerate(paths_list)]...)
#             M2_2 = sort(sort.(M2_2))
#             M2_2_counter = counter(M2_2)
#             M2_2 = []
#             for (key, val) in M2_2_counter
#                 if mod(val, 2)≠0
#                     push!(M2_2, key)
#                 end
#             end
            

#             push!(paths_list2, sort(sort.(M2_1)))
#             push!(paths_list2, sort(sort.(M2_2)))
#         end
#         paths_list2 = collect(Set(paths_list2))
#         dists_list2 = [sum([g[e[1], e[2]] for e in M]) for M in paths_list2]
#         ind = sortperm(dists_list2)
#         dists_list2 = dists_list2[ind]
#         paths_list2 = paths_list2[ind]
        
#         M1 = mwm(g, syndrome)
#         w_M1 = sum([g[e[1], e[2]] for e in M1])
#         if !(w_M1 ≈ dists_list2[1])
#             println(g)
#             println(syndrome)
#         end
#         @test w_M1 ≈ dists_list2[1]
#         @test paths_list2[1] == M1
        
#         wc, c = minimum_weight_cycle(g)

#         if dists_list2[1] + wc < dists_list2[2]
#             println("mwc needed")
#             M2 = vcat(M1, [[c[j], c[j+1]] for j in 1 : length(c)-1])
#             M2 = sort(sort.(M2))
#         else
#             M2 = paths_list2[2]
#         end
#         w_M2 = sum([g[e[1], e[2]] for e in M2])
#         if !(w_M2_true ≈ w_M2)
#             println(g)
#             println(syndrome)
#         end
#         @test M2_true == M2
#         @test w_M2_true ≈ w_M2

#         # Add the tests for mwms with K=2
#         XK = mwms(g, syndrome, 2)
#         @test XK[1] == M1
#         @test XK[2] == M2_true

#         # Test the new data structure for g
#         g2 = get_graph_from_matrix(g)
#         XK2 = mwms(g2, syndrome, 2)
#         M1_2, M2_3 = XK2[1], XK2[2]
#         M1_2 = [g2[e][1] for e in M1_2]
#         M2_3 = [g2[e][1] for e in M2_3]
#         @test sort.(sort.(M1_2)) == M1
#         @test sort.(sort.(M2_3)) == M2_true
#     end
# end

# ## Test secondmwm for larger K for arbitrary graphs with two syndromes

# test_cases = [[2, 5], [1, 3], [4, 8], [7, 6]]
# num_samples = 1000
# N = 8 
# K = 10
# for _ in 1:num_samples
#     g = rand(N, N) * N/2
#     g = g + transpose(g)
#     g = g - diagm(diag(g))

#     for case in test_cases
#         source, target = case[1], case[2]
#         res1 = shortest_paths(g, source, target, K)
#         dists, previous_shortest_paths = res1[1], res1[2]

#         wc, c = minimum_weight_cycle(g)

#         for i in 1 : length(dists)

#         end

#         if res1[1][2] > wc + res1[1][1] 
#             # the second mwm is the shortest path with the mwc
#             M2_v2 = [res1[2][1], [[c[j], c[j+1]] for j in 1 : length(c)-1]...]
            
#         else
#             # the second mwm is the 2nd shortest path
#             M2_v2 = res1[2][2]
#             M2_v2 = [[M2_v2[j], M2_v2[j+1]] for j in 1 : length(M2_v2)-1]
#         end
#         M2_v2 = sort(sort.(M2_v2))
#         w_M2_v2 = sum([g[e[1], e[2]] for e in M2_v2])

#         syndrome = zeros(Int, size(g, 1))
#         syndrome[source] = syndrome[target] = 1
#         w_M2, M2 = secondmwm(g, syndrome)

#         if !(w_M2 ≈ w_M2_v2)
#             println(g)
#             println(syndrome)
#         end        
#         @test w_M2 ≈ w_M2_v2
#         @test M2 == M2_v2

#         # Test the new data structure for g
#         g2 = get_graph_from_matrix(g)
#         w_M2_2, M2_2 = secondmwm(g2, syndrome)
#         M2_2 = sort(sort.([g2[e][1] for e in M2_2]))
#         @test w_M2 ≈ w_M2_2
#         @test M2 == M2_2
#     end
# end

# ## Test mwms for larger K give nondecreasing weights for arbitrary graphs with multiple syndromes

# num_samples = 100
# Nrange = [8, 10, 12]
# num_mwms = 100
# for N in Nrange
#     println("N=$N")
#     for _ in 1:num_samples 
#         g = rand(N, N) * N/2
#         g = g + transpose(g)
#         g = g - diagm(diag(g))

#         syndrome = round.(Int, rand(size(g, 1)-1))
#         if mod(sum(syndrome), 2) == 1
#             push!(syndrome, 1)
#         else
#             push!(syndrome, 0)
#         end
        
#         syndrome == zeros(Int, length(syndrome)) && continue
            
#         XK = mwms(g, syndrome, num_mwms)
#         weights = [sum([g[ee[1], ee[2]] for ee in M]) for M in XK] 
#         diff_weights = [weights[i+1] - weights[i] for i in 1 : length(weights)-1]
#         @test all(diff_weights .>= 0)
#         @test length(Set(XK)) == length(XK) # No repeated matching

#         # Test the new data structure for g
#         g2 = get_graph_from_matrix(g)
#         XK2 = mwms(g2, syndrome, num_mwms)
#         weights2 = [sum([g2[ee][2] for ee in M]) for M in XK2] 
#         diff_weights = [weights[i+1] - weights[i] for i in 1 : length(weights)-1]
#         @test all(diff_weights .>= 0)
#         @test length(Set(XK2)) == length(XK2) # No repeated matching
#     end
# end

# ## Test mwms for larger K for arbitrary graphs with multiple syndromes

# function find_K_combinations(dists_list)
#     """
#     Given a dists_list with n rows, each row has K element.
#     We need to pick one for each row to make a combination
#     We want to find the first K min weight combination
#     There are in total K^n combinations
#     We need a priority queue
#     """
    
#     n = length(dists_list)
#     K = length(dists_list[1])

#     XK = [ones(Int, n)]
#     X = PriorityQueue()
#     X[ones(Int, n)] = sum([dists[1] for dists in dists_list])
#     for k in 2 : K
#         cp = dequeue!(X)
#         for j in 1 : n
#             cp2 = deepcopy(cp)
#             if cp2[j] == K
#                 continue
#             else
#                 cp2[j] += 1
#                 if cp2 ∉ XK
#                     w_cp2 = sum([dists[ind] for (dists, ind) in zip(dists_list, cp2)])
#                     X[cp2] = w_cp2
#                 end
#             end
#         end
#         push!(XK, peek(X)[1])
#     end

#     combs = [[dists_list[i][x[i]] for i in 1:length(x)] for x in XK]
#     weights = [sum(i) for i in combs]
    
#     # return weights, combs
#     return weights, XK
# end
# testcase = [[1,2,3,4], [1,2,3,4]]
# @test find_K_combinations(testcase) == ([2, 3, 3, 4], [[1, 1], [2, 1], [1, 2], [2, 2]])
# testcase = [[1,2,3,4], [3,4,5,6]]
# @test find_K_combinations(testcase) == ([4, 5, 5, 6], [[1, 1], [2, 1], [1, 2], [2, 2]])


# Nrange = [8, 10, 12, 14]
# num_samples = 100
# num_mwms = 20

# for N in Nrange
#     println("$N")
#     for _ in 1 : num_samples
#         g = rand(N, N) * N/2
#         g = g + transpose(g)
#         g = g - diagm(diag(g))

#         syndrome = round.(Int, rand(size(g, 1)-1))
#         if mod(sum(syndrome), 2) == 1
#             push!(syndrome, 1)
#         else
#             push!(syndrome, 0)
#         end

#         syndrome == zeros(Int, length(syndrome)) && continue

#         case = findall(syndrome .== 1)
#         subcases = partition_array(case)

#         M_list = []
#         for subcase in subcases
#             dists_list = []
#             paths_list = []
#             for subsubcase in subcase
#                 source, target = subsubcase[1], subsubcase[2]
#                 dists, paths = shortest_paths(g, source, target, num_mwms)
#                 paths = [[[path[j], path[j+1]] for j in 1 : length(path)-1] for path in paths]
#                 push!(dists_list, dists)
#                 push!(paths_list, paths)
#             end

#             _, Ms_ind = find_K_combinations(dists_list)
#             Ms = [vcat([paths_list[i][x[i]] for i in 1:length(x)]...) for x in Ms_ind]
#             # Ms = [sort(sort.(M)) for M in Ms]
#             # Ms2 = []
#             # M_list = vcat(M_list, Ms)
#             for M in Ms
#                 M = sort(sort.(M))
#                 matching_counter = counter(M)
#                 matching = Vector{Vector{Int64}}()
#                 for (key, val) in matching_counter
#                     if mod(val, 2)≠0
#                         push!(matching, key)
#                     end
#                 end
#                 if length(matching) > 0 
#                     matching = sort(sort.(matching))
#                     w_matching = sum([g[e[1], e[2]] for e in matching])
#                     if matching ∉ M_list
#                         # println(matching)
#                         push!(M_list, matching)
#                     end
#                 end
#             end
#         end
#         weight_list = [sum([g[e[1], e[2]] for e in M]) for M in M_list]
#         ind = sortperm(weight_list)
#         weight_list = weight_list[ind]
#         M_list = M_list[ind]

#         if length(M_list) > num_mwms
#             M_list = M_list[1:num_mwms]
#             weight_list = weight_list[1:num_mwms]
#         end

#         cycles = mwms(g, zeros(Int, size(g, 1)), num_mwms)
#         w_cycles = [sum([g[e[1], e[2]] for e in M]) for M in cycles]

#         new_M_list = []
#         new_weights_list = []
#         for (M, weight) in zip(M_list, weight_list)
#             for (cycle, w) in zip(cycles, w_cycles)
#                 if w + weight > weight_list[end]
#                     break
#                 else
#                     matchings = [M, cycle]
#                     matching = vcat(matchings...)
#                     matching = sort(sort.(matching))
#                     matching_counter = counter(matching)
#                     matching = Vector{Vector{Int64}}()
#                     for (key, val) in matching_counter
#                         if mod(val, 2)≠0
#                             push!(matching, key)
#                         end
#                     end
#                     if length(matching) > 0 
#                         matching = sort(sort.(matching))
#                         w_matching = sum([g[e[1], e[2]] for e in matching])
#                         if matching ∉ M_list && matching ∉ new_M_list && w_matching < weight_list[end]
#                             # println(matching)
#                             push!(new_M_list, matching)
#                             push!(new_weights_list, w_matching)
#                         end
#                     end
#                 end
#             end
#         end

#         weight_list = vcat(weight_list, new_weights_list)
#         M_list = vcat(M_list, new_M_list)
#         ind2 = sortperm(weight_list)
#         weight_list = weight_list[ind2]
#         M_list = M_list[ind2]
#         if length(weight_list) > num_mwms
#             weight_list, M_list = weight_list[1:num_mwms], M_list[1:num_mwms]
#         end


#         M_list_true = mwms(g, syndrome, num_mwms)

#         weight_list_true = [sum([g[e[1], e[2]] for e in M]) for M in M_list_true]
#         # println(length.([weight_list_true, weight_list]))
#         if !(weight_list_true ≈ weight_list)
#             println(g)
#             println(syndrome)
#         end
#         @test weight_list_true ≈ weight_list
#         @test M_list_true == M_list

#         # Test the new data structure for g
#         g2 = get_graph_from_matrix(g)
#         M_list_2 = mwms(g2, syndrome, num_mwms)
#         weight_list_2 = [sum([g2[e][2] for e in M2]) for M2 in M_list_2]
#         M_list_2 = [sort.(sort.([g2[e][1] for e in M2])) for M2 in M_list_2]
#         @test weight_list_2 ≈ weight_list_true
#         @test M_list_2 == M_list_true
#     end
# end
