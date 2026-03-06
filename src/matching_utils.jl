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

"""
    function mwpm(g::Matrix, syndrome::Vector{Int64})

Given the adjacency matrix of a weighted graph g, and a syndrome, 
return the minimum-weight-perfect-matching (MWPM).

Here g is an nxn square matrix and syndrome is an n-component vector that 
sums to an even number. The syndrome indicates a set of even number of 
vertices, and the return mwpm contains a set of (weighted) edges that 
matche the highlighted vertices in pair with minimized total weights. 
"""
function mwpm(g::Matrix, syndrome::Vector{Int64}; ϵ = 1e-16)

    if !(g ≈ transpose(g))
        error("The adjacency matrix of the graph has to be symetric!")
    end

    if !(length(syndrome) == size(g, 1))
        error("The length of the syndrome has to be the same as the number of rows in the adjacency matrix!")
    end

    if !(mod(sum(syndrome), 2) == 0)
        error("The sum of the syndrome has to be even!")
    end
    
    
    fault_id = 0 # 0-base for python
    vertices_edge_mapping = Dict{Int64, Tuple{Int64, Int64}}()
    G = pymatching.Matching()
    for i in range(1, size(g, 1))
        for j in range(i+1, size(g, 2))
            if abs(g[i, j]) > ϵ
                G.add_edge(i-1, j-1, weight=g[i, j], fault_ids=@py {fault_id})
                vertices_edge_mapping[fault_id+1] = (i, j) # 1-base for julia
                fault_id += 1
            end
        end
    end

    selected_edges = pyconvert(Vector{Int}, G.decode(syndrome))

    # return selected_edges, vertices_edge_mapping
    path_dict = filter(p->selected_edges[p.first]==1, vertices_edge_mapping)
    return values(path_dict)
end


# The above the mwpm used for the closest point paper. We won't modify it

###############################################################################################
################ Functions below are for finding K mwms for a given graph ####################
###############################################################################################


function decoding_graph(ηs::Vector, stabilizers::Dict{Int64, Vector{Int64}}, η::Vector; ϵ=1e-10)

    x = ηs/√π ./η 

    # get the weights for each qubit, which are the edges of the graph
    closest_integers = closest_integer.(x)
    second_closest_integers = second_closest_integer.(x)

    edge_weight_list = ((second_closest_integers.-x).^2 - (closest_integers.-x).^2) .* η.^2

    # get edge_to_mode_dict
    num_real_vertices, num_vertices = length(stabilizers), length(stabilizers)
    keys_stabilizers = collect(keys(stabilizers))
    values_stabilizers = collect(values(stabilizers))
    edge_to_mode_dict = Dict()
    for qubit in 1 : length(x)
        qubit_in_stab = findall(qubit .∈ values_stabilizers) # determine which stabilizers the qubit is in
        if length(qubit_in_stab) > 2
            error("Cannot decode code where a single fault can lead more than 2 errors.")
        elseif length(qubit_in_stab) == 2
            vertex_1, vertex_2 = keys_stabilizers[qubit_in_stab[1]], keys_stabilizers[qubit_in_stab[2]]
        elseif length(qubit_in_stab) == 1
            num_vertices += 1                        
            vertex_1, vertex_2 = keys_stabilizers[qubit_in_stab[1]], num_vertices
        elseif length(qubit_in_stab) == 0 # ignore if the qubit is not in any stabilizers
            continue
        end

        merge!(edge_to_mode_dict, Dict(Set([vertex_1, vertex_2]) => (qubit, edge_weight_list[qubit])))
    end

    # highlight the unhappy stabilizers/vertices
    highlighted_vertices = zeros(Int, num_vertices)
    for (index, stabilizer) in stabilizers
        if mod(sum(closest_integers[stabilizer]), 2) == 1
            highlighted_vertices[index] = 1
        end
    end

    highlighted_vertices[num_real_vertices+1] = mod(sum(highlighted_vertices), 2)

    rows = Vector{Int64}()
    cols = Vector{Int64}()
    vals = Vector{Float64}()
    for (key, (_, weight)) in edge_to_mode_dict
        vertex_1, vertex_2 = collect(key)[1], collect(key)[2]        
        push!(rows, vertex_1)
        push!(cols, vertex_2)
        push!(vals, weight)
        push!(rows, vertex_2)
        push!(cols, vertex_1)
        push!(vals, weight)
    end
    for i in num_real_vertices+1 : num_vertices - 1
        push!(rows, i)
        push!(cols, i+1)
        push!(vals, ϵ)
        push!(rows, i+1)
        push!(cols, i)
        push!(vals, ϵ)        
    end    
    g = sparse(rows, cols, vals, num_vertices, num_vertices)        
    
    fault_id = 0 # 0-base for python
    vertices_edge_mapping = Dict{Int64, Tuple{Int64, Int64}}()
    G = pymatching.Matching()

    for (i, j, val) in zip(findnz(g)...)
        if i < j
            G.add_edge(i-1, j-1, weight=val, fault_ids= fault_id)
            vertices_edge_mapping[fault_id+1] = (i, j) # 1-base for julia
            fault_id += 1
        end
    end
    return g, highlighted_vertices, edge_to_mode_dict, closest_integers, second_closest_integers, G, vertices_edge_mapping
end

decoding_graph(ηs::Vector, stabilizers::Dict{Int64, Vector{Int64}}) = decoding_graph(ηs::Vector, stabilizers::Dict{Int64, Vector{Int64}}, ones(length(ηs)))        


function shortest_path(
    g::SparseMatrixCSC{T},
    source::U,
    target::U
) where {T<:Real} where {U<:Integer}

    nvg = size(g, 1)

    if source > size(g, 1)
        error("`source` cannot be larger than the number of vertices in the graph.")
    end
    if target > size(g, 1)
        error("`target` cannot be larger than the number of vertices in the graph.")
    end

    if source == target
        return 0, [source]
    end


    dists = fill(typemax(T), nvg) # Distance of the node to the source
    parents = zeros(U, nvg) # Parent of the nodes
    visited = zeros(Bool, nvg) # If the nodes have been visited

    H = PriorityQueue{U,T}() # Priority queue of the visited nodes, sorted by their distances

    dists[source] = zero(T) # The distance between the source and itself is zero
    visited[source] = true # The source node has been visited
    H[source] = zero(T) # Set the distance of source in the queue to zero

    while !isempty(H)
        u = dequeue!(H) # Pop the node with min distance

        # Because the shortest path between u and the source has been found
        # if we have reached the target node, exit
        if u == target 
            break
        end


        d = dists[u] # The dist of u to the source # Cannot be typemax if `u` is in the queue
        outneighbors2 = findall(g[u,:] .> 0) # The outgoing neighbors for u [note that we also considered directed graphs]
        for v in outneighbors2
            alt = d + g[u, v] # The distance of v to the source

            if !visited[v] 
                # If v has not been visited, then 
                visited[v] = true # mark it as visisted
                dists[v] = alt # Update its distance
                parents[v] = u # Update its parent
                H[v] = alt # Add it into the queue
            elseif alt < dists[v]
                # If v has been visited, but previous distance is larger than the current one 
                dists[v] = alt # Update its distance
                parents[v] = u # Update its parent
                H[v] = alt # Update the distance of v in the queue
            elseif alt >= dists[v] # If v has been visited, and the previous distance is smaller, do nothing
            end
        end
    end

    parents[source] = 0

    function spath(target, parents, source)
        return if target == 0
            nothing
        elseif target == source
            target
        else
            [spath(parents[target], parents, source) target]
        end
    end    
    path = spath(target, parents, source)

    if nothing in path
        return Inf, []
    else
        return dists[target], vec(path)
    end
end

function minimum_weight_cycle(g::SparseMatrixCSC)
    weight = 1e10
    edges = findall(triu(g, 1).!=0) # exclude the main diagonal
    cycle = nothing
    for edge in edges
        weight_edge = g[edge[1], edge[2]]
        g[edge[1], edge[2]] = 0
        g[edge[2], edge[1]] = 0
        dist, path = shortest_path(g, edge[1], edge[2])

        if weight_edge + dist < weight
            weight = weight_edge + dist
            cycle = path
            push!(cycle, cycle[1])
        end
        g[edge[1], edge[2]] = g[edge[2], edge[1]] = weight_edge
    end

    return weight, cycle
end



function mwm(g, G, syndrome, vertices_edge_mapping; ϵ=1e-16, option=2)
#         myid()==2 && println("i am here")
    if syndrome == zeros(Int, length(syndrome))
        c = minimum_weight_cycle(g)[2] # mwc in the format of [1,2,3,1]
        if c == nothing
            M1 = nothing
        else
            M1 = [[c[j], c[j+1]] for j in 1 : length(c)-1] # mwc in the format of [[1, 2], [2, 3], [3, 1]]
        end

        if M1 == nothing
            return nothing
        else
            return sort(sort.(M1))
        end
    end

    M1 = []
    try
        if option == 1 
            error("no option = 1")
        elseif option == 2
            selected_edges = pyconvert(Vector{Int}, G.decode(syndrome))
            path_dict = filter(p->selected_edges[p.first]==1, vertices_edge_mapping)
        end
        M1 = values(path_dict)
        M1 = sort(sort.(collect.(M1))) 
            
        syndrome_2 = deepcopy(syndrome)
#             println([M1, syndrome_2])
        for e in M1
            syndrome_2[e] = 1 .- syndrome_2[e]
        end
            
        if sum(syndrome_2) == 0
            return M1
        else
            return nothing
        end                
    catch
        # If M1 cannot be defined, it means there is no way to match the 
        # highlighted vertices. So we return nothing
        # For example, if the graph G has two disjoint parts G1, and G2, 
        # there is no M1 if the two highlighted vertices are in G1 and G2
        # respectively.
        return nothing
    end  
end

function mwms(g, syndrome::Vector{U}, G, vertices_edge_mapping, K::U; ϵ::T = 1e-16, Δ::T=1e4) where {T<:Real} where {U<:Integer}
    if K<=0
        error("`K` has to be an positive intger.")
    end

    M1 = mwm(g, G, syndrome, vertices_edge_mapping; ϵ = ϵ, option = 2)
    M1 == nothing && return []

    XK = [M1]
    w_M1 = sum([g[ee[1], ee[2]] for ee in M1])

    X = PriorityQueue{Tuple, Real}()
    X[(M1, [], syndrome, [])] = w_M1
    for k in 2 : K

        MK, gp0, syndromep, edgespp = dequeue!(X) # Pop the matching with min weight
        M1p = setdiff(MK, edgespp)
        tt0 = 0
        for j in 1 : length(M1p)+1
#                 println([k, j])
            syndrome_j = deepcopy(syndromep)

            gp = deepcopy(gp0)
            deleted_values = []
            fault_ids = []
            for e in gp
                push!(deleted_values, g[e[1], e[2]])
                g[e[1], e[2]] = g[e[2], e[1]] = 0

                push!(fault_ids, G.get_edge_data(e[1]-1, e[2]-1)["fault_ids"])
                G.add_edge(e[1]-1, e[2]-1, weight=Δ, fault_ids = -1, merge_strategy="replace")
            end

            for l in 1 : min(j, length(M1p))
                e = M1p[l]
                @assert g[e[1], e[2]] > 0
                push!(deleted_values, g[e[1], e[2]])                
                g[e[1], e[2]] = g[e[2], e[1]] = 0
                push!(gp, e)
                
                @assert pyconvert(Float64, G.get_edge_data(e[1]-1, e[2]-1)["weight"]) > 0
                push!(fault_ids, G.get_edge_data(e[1]-1, e[2]-1)["fault_ids"])
                G.add_edge(e[1]-1, e[2]-1, weight=Δ, fault_ids = -1, merge_strategy="replace")
            end

            for l in 1 : j-1
                e = M1p[l]
                syndrome_j[e] = 1 .- syndrome_j[e]
            end

            M_j = mwm(g, G, syndrome_j, vertices_edge_mapping; ϵ = ϵ, option=2)
                
            for (e, val, fault_id) in zip(gp, deleted_values, fault_ids)
                g[e[1], e[2]] = g[e[2], e[1]] = val
                G.add_edge(e[1]-1, e[2]-1, weight=val, fault_ids = fault_id, merge_strategy="replace")
            end

            if M_j !== nothing
                edgespp_updated = deepcopy(edgespp)

                edgespp_updated = vcat(edgespp_updated, M1p[1:j-1])                
                M_j = vcat(M_j, edgespp_updated)                    

                w_M_j = sum([g[ee[1], ee[2]] for ee in M_j])
                M_j = sort(sort.(M_j))

                if M_j ∉ XK 
                    # Because the matchings have fixed ordering, 
                    # If M_j ∈ X, it will be replaced automatically with no change of weight
                    # Hence there is no need to check if M_j ∈ X
#                     X[(M_j, g_j, syndrome_j, edgespp_updated)] = w_M_j
#                         tt += @elapsed 
                    X[(M_j, gp, syndrome_j, edgespp_updated)] = w_M_j
                end
            end
        end
#             push!(tt, tt0)
#             myid() == 2 && println("mwms_v5 tt=$tt")
        length(X) == 0 && return Inf, []

        (Mnext, _, _, _), w_Mnext = peek(X)
        push!(XK, Mnext)

        w_Mnext == Inf && break
    end
#         myid()==2 && println("v5 $(length(tt))")
    
    XK2 = []
    while length(X) > 0
        MK, gp0, syndromep, edgespp = dequeue!(X)
        push!(XK2, MK)
    end
    
    return XK, XK2 #, tt, t_mwc
end

## Update decoding_graph, mwm and mwms
function decoding_graph_v2(ηs::Vector, stabilizers::Dict{Int64, Vector{Int64}}, η::Vector; ϵ=1e-10)

    x = ηs/√π ./η 

    # get the weights for each qubit, which are the edges of the graph
    closest_integers = closest_integer.(x)
    second_closest_integers = second_closest_integer.(x)

    edge_weight_list = ((second_closest_integers.-x).^2 - (closest_integers.-x).^2) .* η.^2

    # get edge_to_mode_dict
    num_real_vertices, num_vertices = length(stabilizers), length(stabilizers)
    keys_stabilizers = collect(keys(stabilizers))
    values_stabilizers = collect(values(stabilizers))
    edge_to_mode_dict = Dict()
    for qubit in 1 : length(x)
        qubit_in_stab = findall(qubit .∈ values_stabilizers) # determine which stabilizers the qubit is in
        if length(qubit_in_stab) > 2
            error("Cannot decode code where a single fault can lead more than 2 errors.")
        elseif length(qubit_in_stab) == 2
            vertex_1, vertex_2 = keys_stabilizers[qubit_in_stab[1]], keys_stabilizers[qubit_in_stab[2]]
        elseif length(qubit_in_stab) == 1
            num_vertices += 1                        
            vertex_1, vertex_2 = keys_stabilizers[qubit_in_stab[1]], num_vertices
        elseif length(qubit_in_stab) == 0 # ignore if the qubit is not in any stabilizers
            continue
        end

        merge!(edge_to_mode_dict, Dict(Set([vertex_1, vertex_2]) => (qubit, edge_weight_list[qubit])))
    end

    # highlight the unhappy stabilizers/vertices
    highlighted_vertices = zeros(Bool, num_vertices)
    for (index, stabilizer) in stabilizers
        if mod(sum(closest_integers[stabilizer]), 2) == 1
            highlighted_vertices[index] = 1
        end
    end

    highlighted_vertices[num_real_vertices+1] = mod(sum(highlighted_vertices), 2)

    rows = Vector{Int64}()
    cols = Vector{Int64}()
    vals = Vector{Float64}()
    for (key, (_, weight)) in edge_to_mode_dict
        vertex_1, vertex_2 = collect(key)[1], collect(key)[2]        
        push!(rows, vertex_1)
        push!(cols, vertex_2)
        push!(vals, weight)
        push!(rows, vertex_2)
        push!(cols, vertex_1)
        push!(vals, weight)
    end
    for i in num_real_vertices+1 : num_vertices - 1
        push!(rows, i)
        push!(cols, i+1)
        push!(vals, ϵ)
        push!(rows, i+1)
        push!(cols, i)
        push!(vals, ϵ)        
    end    
    g = sparse(rows, cols, vals, num_vertices, num_vertices)        

    fault_id = 0 # 0-base for python
    vertices_edge_mapping = Dict{Int64, Vector{Int64}}()
    edge_vertices_mapping = Dict{Vector{Int64}, Int64}()
    G = pymatching.Matching()

    for (i, j, val) in zip(findnz(g)...)
        if i < j
            G.add_edge(i-1, j-1, weight=val, fault_ids= fault_id)
            vertices_edge_mapping[fault_id+1] = [i, j] # 1-base for julia
            edge_vertices_mapping[[i, j]] = fault_id+1 # 1-base for julia
            fault_id += 1
        end
    end
    return g, highlighted_vertices, edge_to_mode_dict, closest_integers, second_closest_integers, G, vertices_edge_mapping, edge_vertices_mapping
end


decoding_graph_v2(ηs::Vector, stabilizers::Dict{Int64, Vector{Int64}}) = decoding_graph_v2(ηs::Vector, stabilizers::Dict{Int64, Vector{Int64}}, ones(length(ηs)))        


function mwm_v2(g::SparseMatrixCSC, G, syndrome::Vector{Bool}, vertices_edge_mapping::Dict{Int64, Vector{Int64}}, edge_vertices_mapping::Dict{Vector{Int64}, Int64})
#         myid()==2 && println("i am here")
    if sum(syndrome) == 0 # syndrome == zeros(Bool, length(syndrome))
        c = minimum_weight_cycle(g)[2] # mwc in the format of [1,2,3,1]
        if c == nothing
            return nothing
        else
#                 M1 = [[c[j], c[j+1]] for j in 1 : length(c)-1] # mwc in the format of [[1, 2], [2, 3], [3, 1]]
#                 M1 = [edge_vertices_mapping[min(c[j], c[j+1]), max(c[j], c[j+1])] for j in 1 : length(c)-1]
            M1 = [edge_vertices_mapping[sort(c[j:j+1])] for j in 1 : length(c)-1]
            return sort(M1)
        end
    end

    M1 = []
    try
        selected_edges = pyconvert(Vector{Bool}, G.decode(syndrome))
        M1 = findall(selected_edges)
        syndrome_2 = deepcopy(syndrome)
        for e in M1
            e2 = vertices_edge_mapping[e]
            syndrome_2[e2] = 1 .- syndrome_2[e2]
        end

        if sum(syndrome_2) == 0
            return sort(M1)
        else
            return nothing
        end                
    catch
        # If M1 cannot be defined, it means there is no way to match the 
        # highlighted vertices. So we return nothing
        # For example, if the graph G has two disjoint parts G1, and G2, 
        # there is no M1 if the two highlighted vertices are in G1 and G2
        # respectively.
        return nothing
    end  
end

function get_next_matching!(g::SparseMatrixCSC, X, XK, G, vertices_edge_mapping::Dict{Int64, Vector{Int64}}, edge_vertices_mapping::Dict{Vector{Int64}, Int64}; Δ::Real=1e2)
    MK, gp0, syndromep, edgespp = dequeue!(X) # Pop the matching with min weight
    M1p = setdiff(MK, edgespp)
    for j in 1 : length(M1p)+1
        syndrome_j = deepcopy(syndromep)

        gp = deepcopy(gp0)
        deleted_values = []
        fault_ids = []
        for e in gp
            e2 = vertices_edge_mapping[e]
            push!(deleted_values, g[e2...])
            g[e2...] = g[reverse(e2)...] = Δ
            fault_id = G.get_edge_data((e2 .- 1)...)["fault_ids"]
            push!(fault_ids, fault_id)
            G.add_edge((e2 .- 1)..., weight=Δ, fault_ids = fault_id, merge_strategy="replace")
        end

        for l in 1 : min(j, length(M1p))
            e = M1p[l]
            e2 = vertices_edge_mapping[e]
            @assert g[e2...] > 0
            push!(deleted_values, g[e2...])
            g[e2...] = g[reverse(e2)...] = Δ
            push!(gp, e)

            @assert pyconvert(Float64, G.get_edge_data((e2 .- 1)...)["weight"]) > 0
            fault_id = G.get_edge_data((e2 .- 1)...)["fault_ids"]
            push!(fault_ids, fault_id)
            G.add_edge((e2 .- 1)..., weight=Δ, fault_ids = fault_id, merge_strategy="replace")
        end

        for l in 1 : j-1
            e = M1p[l]
            e2 = vertices_edge_mapping[e]
            syndrome_j[e2] = 1 .- syndrome_j[e2]                
        end

        M_j = mwm_v2(g, G, syndrome_j, vertices_edge_mapping, edge_vertices_mapping)

        if M_j !== nothing
            w_M_j = sum([g[vertices_edge_mapping[ee]...] for ee in M_j])
            if w_M_j >= Δ
                M_j = nothing
            end
        end                

        for (e, val, fault_id) in zip(gp, deleted_values, fault_ids)
            e2 = vertices_edge_mapping[e]
            g[e2...] = g[reverse(e2)...] = val
            G.add_edge((e2 .- 1)..., weight=val, fault_ids = fault_id, merge_strategy="replace")                
        end

        if M_j !== nothing
            edgespp_updated = deepcopy(edgespp)

            # no need to sort edgespp_updated because M1p is already sorted
            edgespp_updated = vcat(edgespp_updated, M1p[1:j-1])                
            M_j = vcat(M_j, edgespp_updated)                    

            w_M_j = sum([g[vertices_edge_mapping[ee]...] for ee in M_j])
            M_j = sort(M_j)

            if M_j ∉ XK 
                # Because the matchings have fixed ordering, 
                # If M_j ∈ X, it will be replaced automatically with no change of weight
                # Hence there is no need to check if M_j ∈ X
                X[(M_j, gp, syndrome_j, edgespp_updated)] = w_M_j
            end
        end
    end
    return X
end

function mwms_v2(g::SparseMatrixCSC, syndrome::Vector{Bool}, G, vertices_edge_mapping::Dict{Int64, Vector{Int64}}, edge_vertices_mapping::Dict{Vector{Int64}, Int64}, K::Int; Δ::Real=1e2, option::Int=1)
    if K<=0
        error("`K` has to be an positive intger.")
    end
    if option ∉ [1, 2]
        error("`option` can only be 1 or 2.")
    end        
    
    if sum(syndrome) == 0 && K==1
        if option == 1
            return [[]], []
        else
            X = PriorityQueue{Tuple, Real}()
            X[([], [], syndrome, [])] = 0
            return [[]], X
        end
    end     

    M1 = mwm_v2(g, G, syndrome, vertices_edge_mapping, edge_vertices_mapping)
    M1 == nothing && return []

    XK = [M1]
    w_M1 = sum([g[vertices_edge_mapping[ee]...] for ee in M1])

    X = PriorityQueue{Tuple, Real}()
    X[(M1, [], syndrome, [])] = w_M1

    if sum(syndrome) == 0
        pushfirst!(XK, [])
        Kloop = K-1
    else
        Kloop = K
    end    

    for k in 2 : Kloop
        get_next_matching!(g, X, XK, G, vertices_edge_mapping, edge_vertices_mapping; Δ=Δ)
        length(X) == 0 && return Inf, []
        (Mnext, _, _, _), w_Mnext = peek(X)
        push!(XK, Mnext)
        w_Mnext == Inf && break
    end
    if option == 1
        XK2 = []
        while length(X) > 0
            MK, gp0, syndromep, edgespp = dequeue!(X)
            push!(XK2, MK)
        end
        return XK, XK2
    else
        return XK, X
    end
end


## The following functions are mainly for testing the functions above

shortest_path(
    g::Matrix{T},
    source::U,
    target::U
) where {T<:Real} where {U<:Integer} = shortest_path(sparse(g), source, target)

minimum_weight_cycle(g::Matrix{T}) where {T<:Real} = minimum_weight_cycle(sparse(g))

"""
    function next_shortest_path!(
        g::Matrix, 
        source::Int64, 
        target::Int64, 
        previous_shortest_paths::Vector,
        potential_shortest_paths::Vector; 
        directed_graph::Bool = false
    )

Given the adjacency matrix of a weighted graph `g``, and two vertices `source` and `target`, and 
previously found and potential shortest paths, find the next shortest path between the two vertices.

Notes: Copied and modified from https://github.com/JuliaGraphs/Graphs.jl/blob/master/src/shortestpaths/yen.jl

Notes: The graph cannot have negative weighted cycle

Notes: previous_shortest_paths and potential_shortest_paths are modified in place
"""
function next_shortest_path!(
    g::Matrix{T},
    source::U,
    target::U,
    previous_shortest_paths::Vector{Vector{U}},
    potential_shortest_paths::PriorityQueue{Vector{U}, T},
    directed_graph::Bool = false
) where {T<:Real} where {U<:Integer}

    if length(previous_shortest_paths) == 0
        dist, path = shortest_path(g, source, target)
        push!(previous_shortest_paths, path)
        return dist, path
    end

    for j in 1 : length(previous_shortest_paths[end])
        gcopy = deepcopy(g)

        # Spur node is retrieved from the previous k-shortest path, k − 1
        spurnode = previous_shortest_paths[end][j]
        #  The sequence of nodes from the source to the spur node of the previous k-shortest path
        rootpath = previous_shortest_paths[end][1:j]

        # Remove the links of the previous shortest paths which share the same root path
        for ppath in previous_shortest_paths
            if length(ppath) > j && rootpath == ppath[1:j]
                u = ppath[j]
                v = ppath[j + 1]
                gcopy[u, v] = 0
                if directed_graph != true
                    gcopy[v, u] = 0
                end
            end
        end
        
        # Remove node of root path and calculate dist of it
        distrootpath = zero(T)
        for n in 1:(length(rootpath) - 1)
            u = rootpath[n]
            # Evaluate distance of root path
            v = rootpath[n + 1]
            distrootpath += gcopy[u, v]

            if directed_graph == true
                gcopy[u,u:end] .= 0
            else
                gcopy[u, :] .= gcopy[:, u] .= 0
            end
        end
        
        # Calculate the spur path from the spur node to the sink
        dist_spur, sp_spur = shortest_path(gcopy, spurnode, target)
        if length(sp_spur) >= 2
            # Entire path is made up of the root path and spur path
            pathtotal = [rootpath[1:(end - 1)]; sp_spur]
            distpath = distrootpath + dist_spur
            # Add the potential k-shortest path to the heap
            if !haskey(potential_shortest_paths, pathtotal)
                enqueue!(potential_shortest_paths, pathtotal, distpath)
            end
        end
    end

    if isempty(potential_shortest_paths)
        return Inf, []
    else
        path = peek(potential_shortest_paths)[1]
        dist = peek(potential_shortest_paths)[2]
        push!(previous_shortest_paths, dequeue!(potential_shortest_paths))
        return dist, path
    end
end

"""
    function shortest_paths(g::Matrix, source::Int64, target::Int64, K::Int; directed_graph::Bool = false)

Given the adjacency matrix of a weighted graph `g``, and two vertices `source` and `target`, 
find the K shortest weighted paths between the two vertices.

Notes: Copied and modified from https://github.com/JuliaGraphs/Graphs.jl/blob/master/src/shortestpaths/yen.jl

Notes: The graph cannot have negative weighted cycle
"""
function shortest_paths(
    g::Matrix{T},
    source::U,
    target::U,
    K::U,
    directed_graph::Bool = false
) where {T<:Real} where {U<:Integer}
    nvg = size(g, 1)

    if K<=0
        error("`K` has to be an positive intger.")
    end

    if source > size(g, 1)
        error("`source` cannot be larger than the number of vertices in the graph.")
    end
    if target > size(g, 1)
        error("`target` cannot be larger than the number of vertices in the graph.")
    end

    source == target && return [0], [[source]]

    previous_shortest_paths = Vector{Vector{U}}()
    potential_shortest_paths = PriorityQueue{Vector{U}, T}()
    dists = Array{T,1}()
    for k in 1 : K
        dist, _ = next_shortest_path!(g, source, target, previous_shortest_paths, potential_shortest_paths, directed_graph)
        dist == Inf && break
        push!(dists, dist)
    end

    return dists, previous_shortest_paths
end
