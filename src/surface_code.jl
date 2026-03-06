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
    surface_code_Z_stabilizers(d::Int)

Return the dictionary for the Z stabilizers of the surface code with distance d.

Example: 
    >>> surface_code_Z_stabilizers(3)
    >>> Dict{Int64, Vector{Int64}} with 4 entries:
        4 => [6, 9]
        2 => [2, 3, 5, 6]
        3 => [4, 5, 7, 8]
        1 => [1, 4]
"""
function surface_code_Z_stabilizers(d::Int)
    if mod(d, 2) != 1
        error("The distance of a surface code should be odd.")
    end
    Z_dict = Dict{Int64, Vector{Int64}}()
    num_cols = Int((d + 1) / 2)
    num_rows = d - 1
    for row_idx in 1:2:num_rows
        for col_idx in 1:1
            stab_idx = col_idx + (row_idx - 1) * num_cols
            qubit_idx = 2 * col_idx - 2 + (row_idx - 1) * d
            Z_dict = merge!(Z_dict, 
                Dict(stab_idx => [qubit_idx + 1, qubit_idx + d + 1])
            )
        end
        for col_idx in 2:num_cols
            stab_idx = col_idx + (row_idx - 1) * num_cols
            qubit_idx = 2 * col_idx - 2 + (row_idx - 1) * d 
            Z_dict = merge!(Z_dict, 
                Dict(stab_idx => [qubit_idx, qubit_idx + 1, qubit_idx + d, qubit_idx + d + 1])
            )
        end
    end
    for row_idx in 2:2:num_rows
        for col_idx in 1:num_cols-1
            stab_idx = col_idx + (row_idx - 1) * num_cols
            qubit_idx = 2 * col_idx - 1 + (row_idx - 1) * d 
            Z_dict = merge!(Z_dict, 
                Dict(stab_idx => [qubit_idx, qubit_idx + 1, qubit_idx + d, qubit_idx + d + 1])
            )
        end
        for col_idx in num_cols:num_cols
            stab_idx = col_idx + (row_idx - 1) * num_cols
            qubit_idx = 2 * col_idx - 1 + (row_idx - 1) * d 
            Z_dict = merge!(Z_dict, 
                Dict(stab_idx => [qubit_idx, qubit_idx + d])
            )
        end
    end
    return Z_dict
end

"""
    surface_code_X_stabilizers(d::Int)

Return the dictionary for the X stabilizers of the surface code with distance d.

Example: 
    >>> surface_code_X_stabilizers(3)
    >>> Dict{Int64, Vector{Int64}} with 4 entries:
        4 => [5, 6, 8, 9]
        2 => [7, 8]
        3 => [2, 3]
        1 => [1, 2, 4, 5]
"""
function surface_code_X_stabilizers(d::Int)
    if mod(d, 2) != 1
        error("The distance of a surface code should be odd.")
    end
    X_dict = Dict{Int64, Vector{Int64}}()
    num_cols = d - 1
    num_rows = Int((d + 1) / 2)
    for col_idx in 1:2:num_cols
        for row_idx in 1:num_rows-1
            stab_idx = row_idx + (col_idx - 1) * num_rows
            qubit_idx = (row_idx - 1) * 2 * d + col_idx
            X_dict = merge!(X_dict, 
                Dict(stab_idx => [qubit_idx, qubit_idx + 1, qubit_idx + d, qubit_idx + d + 1])
            )
        end
        for row_idx in num_rows:num_rows
            stab_idx = row_idx + (col_idx - 1) * num_rows
            qubit_idx = (row_idx - 1) * 2 * d + col_idx
            X_dict = merge!(X_dict, 
                Dict(stab_idx => [qubit_idx, qubit_idx + 1])
            )
        end
    end
    for col_idx in 2:2:num_cols
        for row_idx in 1:1
            stab_idx = row_idx + (col_idx - 1) * num_rows
            qubit_idx = (row_idx - 1) * 2 * d - d + col_idx
            X_dict = merge!(X_dict, 
                Dict(stab_idx => [qubit_idx + d, qubit_idx + d + 1])
            )
        end
        for row_idx in 2:num_rows
            stab_idx = row_idx + (col_idx - 1) * num_rows
            qubit_idx = (row_idx - 1) * 2 * d - d + col_idx
            X_dict = merge!(X_dict, 
                Dict(stab_idx => [qubit_idx, qubit_idx + 1, qubit_idx + d, qubit_idx + d + 1])
            )
        end
    end
    return X_dict
end

"""
    surface_code_stabilizers(d::Int)

Return the dictionary for the stabilizers of the surface code with distance d.

Example: 
    >>> surface_code_stabilizers(3)
    >>> Dict{Int64, Vector{Int64}} with 8 entries:
      5 => [2, 8]
      4 => [9, 11, 15, 17]
      6 => [4, 6, 10, 12]
      7 => [8, 10, 14, 16]
      2 => [13, 15]
      8 => [12, 18]
      3 => [3, 5]
      1 => [1, 3, 7, 9]
"""
function surface_code_stabilizers(d::Int)
    surface_code_x_stabilizers = surface_code_X_stabilizers(d)
    surface_code_z_stabilizers = surface_code_Z_stabilizers(d)
    
    dict = Dict{Int64, Vector{Int64}}()
    
    for (key, val) in sort(surface_code_x_stabilizers)
        dict[key] = [2item-1 for item in val]
    end

    for (key, val) in sort(surface_code_z_stabilizers)
        dict[length(surface_code_x_stabilizers) + key] = [2item for item in val]
    end    
    
    return dict
end

"""
    surface_code_Z_logicals(d::Int)

Return the Z logical operator for the surface code with distance d

Example: 
    >>> surface_code_Z_logicals(3)
    >>> Dict{Int64, Vector{Int64}} with 1 entry:
        1 => [1, 2, 3]
"""
surface_code_Z_logicals(d::Int) = Dict(1 => collect(1:d))

"""
    surface_code_X_logicals(d::Int)

Return the X logical operator for the surface code with distance d

Example: 
    >>> surface_code_X_logicals(3)
    >>> Dict{Int64, Vector{Int64}} with 1 entry:
        1 => [1, 4, 7]
"""
surface_code_X_logicals(d::Int) = Dict(1 => collect(1:d:d^2))

"""
    surface_code_Mq(d::Int)

Return the generator for the surface code of distance d in the q subspace
"""
function surface_code_Mq(d::Int)
    Mq = 2 * diagm(ones(Int64, d^2))
    X_dict = surface_code_X_stabilizers(d)
    for (_, value_list) in X_dict
        min_value = minimum(value_list)
        for value in value_list
            Mq[min_value, value] = 1
        end
    end
    Mq = Mq/√2 # To make sure det(Mq) = √2
    return Mq
end

"""
    surface_code_Mp(d::Int)

Return the generator for the surface code of distance d in the p subspace
"""
function surface_code_Mp(d::Int)
    Mp = 2 * diagm(ones(Int64, d^2))
    Z_dict = surface_code_Z_stabilizers(d)
    for (_, value_list) in Z_dict
        min_value = minimum(value_list)
        for value in value_list
            Mp[min_value, value] = 1
        end
    end
    Mp = Mp/√2 # To make sure det(Mp) = √2
    return Mp
end

"""
    surface_code_M(d::Int)

Return the generator for the surface code of distance d
"""
function surface_code_M(d::Int)
    Mq = surface_code_Mq(d)
    Mp = surface_code_Mp(d)
    M = BlockDiagonal([Mq, Mp])
    T = basis_transformation(size(Mq)[1])
    return transpose(T) * M * T
end


"""
    decode_surface_rectangular_GKP_code(ξ::Vector, 
    surface_code_x_stabilizers::Dict{String, Dict{Int64, Vector{Int64}}},
    surface_code_z_stabilizers::Dict{String, Dict{Int64, Vector{Int64}}},
    η::Vector
    )

Decode the surface-GKP code for an input vector ξ. 

    Inputs:
        ξ: The input vector
        surface_code_x_stabilizers: The X stabilizers for the stabilizer GKP code
        surface_code_z_stabilizers: The Z stabilizers for the stabilizer GKP code
        η: The square root of aspect ratio for each lattice dimension
    
    Outputs:
        y: The closest point for ξ in the symplectic dual lattice generated by √(2π) * Mperp
        final_list: The logical operator after the attempt of correcting the error
    
    Notes: 
        There are a few assumptions for the surface-GKP code
            1. The inner GKP code has to be one mode rectangular GKP code
            2. The stabilizers should be chosen such that a given shift 
                error can induce at most two syndrome errors
"""
function decode_surface_rectangular_GKP_code(
    ξ::Vector, 
    surface_code_x_stabilizers::Dict{Int64, Vector{Int64}},
    surface_code_z_stabilizers::Dict{Int64, Vector{Int64}},
    η::Vector
)
    d = √(length(ξ)/2)
    if !isinteger(d)
        error("The length of ξ has to be 2d^2 for an integer d")
    end
    
    ξ_q = ξ[1:2:end]
    ξ_p = ξ[2:2:end]
    η_q = η[1:2:end]
    η_p = η[2:2:end]
    

    y_q, final_list_q = decode_concatenated_rectangular_GKP_code(ξ_q, surface_code_z_stabilizers, η_q)
    y_p, final_list_p = decode_concatenated_rectangular_GKP_code(ξ_p, surface_code_x_stabilizers, η_p)

    y = zeros(length(ξ))
    y[1:2:end], y[2:2:end] = y_q, y_p

    final_list = zeros(length(ξ))
    final_list[1:2:end], final_list[2:2:end] = final_list_q, final_list_p

    return y, final_list
end


"""
    decode_surface_code(ξ::Vector, 
    surface_code_x_stabilizers::Dict{String, Dict{Int64, Vector{Int64}}},
    surface_code_z_stabilizers::Dict{String, Dict{Int64, Vector{Int64}}}
    )

Decode the surface square GKP code for an input vector ξ.
"""
decode_surface_code(
    ξ::Vector, 
    surface_code_x_stabilizers::Dict{Int64, Vector{Int64}},
    surface_code_z_stabilizers::Dict{Int64, Vector{Int64}}
) = decode_surface_rectangular_GKP_code(ξ, 
    surface_code_x_stabilizers,
    surface_code_z_stabilizers,
    ones(length(ξ))
)



# An aux function for bsv_surface_code(ηs::Vector, σ::Float64; Nv = 5)
"""
    embed_weights(ws::Vector{Float64}, whichtype::String)

Embed the surface-GKP code with given weights into the unrotated surface code,
and return the corresponding weights in the latter code.

Notes: ws is the weights on the surface-GKP code and whichtype = `x` or `z` 
for X and Z type stabilizers respectively. 
"""
function embed_weights(ws::Vector, whichtype::String)
    d = √length(ws)
    if !isinteger(d)
        error("The length of ws has to be d^2 for an integer d")
    end
    d = Int(d)

    if whichtype == "x"
        GKPind = vcat([((Int((d-1)/2) * (2d-1)+1) : -(d-1) : Int((d+1)/2) ) .+ (i-1) * d for i in 1:d]...)    
    elseif whichtype == "z"
        GKPind = vcat([((Int((d+1)/2)) : d : (Int((d+1)/2)) + d*(d-1)) .+ (i-1) * (d-1) for i in 1:d]...)
    else
        error("`whichtype` can only be `x` or `z`.")
    end

    ws_unrotated = ones(d^2+(d-1)^2)
    ws_unrotated[GKPind] .= ws

    ind0 = []
    # Set the left-top corner data qubits to zero
    for i in 1 : Int((d-1)/2-1)
        ind0 = vcat(ind0, (i-1)*d+1-(i-1)*(d-1) : (d-1) : (i-1)*d+1+(i-1)*(d-1))
    end

    # Set the right-bottom corner data qubits to zero
    for i in 1 : Int((d-1)/2-1)
        ind0 = vcat(ind0, d^2+(d-1)^2-((i-1)*d)-(i-1)*(d-1) : (d-1) : d^2+(d-1)^2-((i-1)*d)+(i-1)*(d-1))
    end
    
    # Set the left-bottom corner data qubits to zero
    for i in 1 : Int((d-1)/2)
        ind0 = vcat(ind0, (i-1)*(d-1)+d-(i-1)*d : d : (i-1)*(d-1)+d+(i-1)*d)
    end
    
    # Set the right-top corner data qubits to zero
    for i in 1 : Int((d-1)/2)
        ind0 = vcat(ind0, d^2+(d-1)^2-d+1-((i-1)*(d-1))-(i-1)*d : d : d^2+(d-1)^2-d+1-((i-1)*(d-1))+(i-1)*d)
    end
    ws_unrotated[ind0] .= 0
    
    return ws_unrotated
end

# An aux function for bsv_surface_code(ηs::Vector, σ::Float64; Nv::Int64 = 5)
"""
    get_prob(ηs::Vector, whichtype::String, σ::Float64; Nv::Int64=5)

Return the log10 of the coset probability given the candidate error ηs of the following

1/√(2πσ^2) * ∑_{g∈G}∏_{j=1}^N∑_{ni=-n}^n exp(-(ηs - √π g_i - 2√π n_i)^2 / (2σ^2))

Notes: wwhichtype = `x` or `z` for X and Z type stabilizers respectively.
σ is the variance of the Gaussian, and n is the number of terms kept in the
summation.

Notes: When g is the identity in G, then the corresponding term reads 
    (∑_{ni=-n}^n exp(-(ηs - 2√π n_i)^2 / (2σ^2)))^N
"""
function get_prob(ηs::Vector, whichtype::String, σ::Float64; Nv::Int64=5)
    ws = []
    w2s = []
    for i in 1 : length(ηs)
        cps = closest_points_Zn([(ηs[i] - √π)/2√π], Nv)
        cp2s = closest_points_Zn([(ηs[i])/2√π], Nv)

        w = sum([exp(-(ηs[i] - √π - cp[1] * 2√π)^2/(2σ^2)) for cp in cps])
        w2 = sum([exp(-(ηs[i] - cp2[1] * 2√π)^2/(2σ^2)) for cp2 in cp2s])

        push!(ws, w)
        push!(w2s, w2)
    end
    
    w3s = ws ./ w2s
    w3s_unrotated = embed_weights(w3s, whichtype)
    Zw3 = bsv(w3s_unrotated) + log10(prod(w2s)) - 1/2 * log10(2π * σ^2)
    
    return Zw3
end

"""
    bsv_surface_code(ηs::Vector, σ::Float64; Nv = 5, subspace="x")

Decode the (rotated) surface code given a candidate error ηs.

Args: 
    ηs: The candidate error
    σ: The noise strength
    Nv: The number of closest points
    subspace: The subspace in which we decode

Return:
    The counter displacement
"""
function bsv_surface_code(ηs::Vector, σ::Float64; Nv = 5, subspace="x")
    if subspace ∈ ["x", "z"]
        d = √(length(ηs))
        correct_length="d^2"
    else
        d = √(length(ηs)/2)
        correct_length="2d^2"
    end

    if !isinteger(d)
        error("The length of ηs has to be $(correct_length) for an integer d")
    end
    d = Int(d)
    
    function bsv_subspace(η, σ, logical, subspace)
        
        # For the q subspace, get the prob with no error and with logical error
        prob_I = get_prob(η, subspace, σ; Nv=Nv)
        prob_X = get_prob(η - logical, subspace, σ; Nv=Nv)

        # Determine if there is logical error
        prob_I < prob_X ? lx = 1 : lx = 0

        if lx == 0
            lstar = zeros(d^2)
        else
            lstar = logical
        end
        return lstar
    end
    
    # Define the X and Z logicals
    x_logical = surface_code_X_logicals(d)[1]
    z_logical = surface_code_Z_logicals(d)[1]

    X = zeros(2d^2)
    X[2x_logical .- 1] .= 1/√2 * √(2π)
    Z = zeros(2d^2)
    Z[2z_logical] .= 1/√2 * √(2π)        
    
    if subspace == "x"
        ηs_q = ηs
        lstar = bsv_subspace(ηs_q, σ, X[1:2:end], subspace)
        return -(ηs_q - lstar)
    elseif subspace == "z"
        ηs_p = ηs
        lstar = bsv_subspace(ηs_p, σ, Z[2:2:end], subspace)
        return -(ηs_p - lstar)
    else
        ηs_q = ηs[1:2:end]
        ηs_p = ηs[2:2:end]
        lstar_q = bsv_subspace(ηs_q, σ, X[1:2:end], "x")
        lstar_p = bsv_subspace(ηs_p, σ, Z[2:2:end], "z")
        rec = ηs
        rec[1:2:end] = rec[1:2:end] - lstar_q
        rec[2:2:end] = rec[2:2:end] - lstar_p
        return -rec
    end
end

"""
    get_coords_surf_hex(dx::Int, dz::Int; shift=0.001)

Return the coordinates for the tensor network for the (rotated) surface code at distance d. The shift is for shifting the coordinates for the tensors of the stabilizers a bit so that the tensor network can be contracted properly when using the SweepContractor.jl package.
"""
function get_coords_surf_hex(dx::Int, dz::Int; shift=0.001)
    if mod(dx, 2) != 1 || dx < 3 || mod(dz, 2) != 1 || dz < 3
        error("The distances of a surface code should be an odd number bigger or equal to 3.")
    end

    coords_data = []
    for i = 0 : dx-1
        for j = 0 : dz-1
            push!(coords_data, (j, -i))
        end
    end

    coords_stab = [[0.0, 0.0] for _ in 1 : dx*dz-1]
    types_stab = ["" for _ in 1 : dx*dz-1]

    # For Z stabs        
    num_cols = Int((dz + 1) / 2)
    num_rows = dx - 1
    num_Z = 0
    for row_idx in 1:2:num_rows
        for col_idx in 1:1
            stab_idx = col_idx + (row_idx - 1) * num_cols
            coords_stab[stab_idx] = [col_idx-3/2 + 1/2 + shift, -(2row_idx-1)/2]
            types_stab[stab_idx] = "Z"                
            num_Z += 1
        end
        for col_idx in 2:num_cols
            stab_idx = col_idx + (row_idx - 1) * num_cols
            coords_stab[stab_idx] = [2col_idx-5/2, -(2row_idx-1)/2]
            types_stab[stab_idx] = "Z"                
            num_Z += 1
        end
    end
    for row_idx in 2:2:num_rows
        for col_idx in 1:num_cols-1
            stab_idx = col_idx + (row_idx - 1) * num_cols
            coords_stab[stab_idx] = [2(col_idx-1)+1/2, -(2row_idx-1)/2]               
            types_stab[stab_idx] = "Z"                
            num_Z += 1
        end
        for col_idx in num_cols:num_cols
            stab_idx = col_idx + (row_idx - 1) * num_cols
            coords_stab[stab_idx] = [2(col_idx-1)+1/2 - 1/2 - shift, -(2row_idx-1)/2]
            types_stab[stab_idx] = "Z"                
            num_Z += 1
        end
    end        

    # For X stabs
    num_cols = dz - 1
    num_rows = Int((dx + 1) / 2)

    for col_idx in 1:2:num_cols
        for row_idx in 1:num_rows-1
            stab_idx = row_idx + (col_idx - 1) * num_rows
            coords_stab[num_Z+stab_idx] = [(col_idx-1)+1/2, -2(row_idx-1)-1/2]
            types_stab[num_Z+stab_idx] = "X"                
        end
        for row_idx in num_rows:num_rows
            stab_idx = row_idx + (col_idx - 1) * num_rows
            coords_stab[num_Z+stab_idx] = [col_idx-1/2, -2(row_idx-1)-1/2+1/2+shift]
            types_stab[num_Z+stab_idx] = "X"                
        end
    end
    for col_idx in 2:2:num_cols
        for row_idx in 1:1
            stab_idx = row_idx + (col_idx - 1) * num_rows
            coords_stab[num_Z+stab_idx] = [col_idx-1/2, 1/2 - 1/2 - shift]
            types_stab[num_Z+stab_idx] = "X"                
        end
        for row_idx in 2:num_rows
            stab_idx = row_idx + (col_idx - 1) * num_rows
            coords_stab[num_Z+stab_idx] = [(col_idx-1)+1/2, -2(row_idx-1)+1/2]
            types_stab[num_Z+stab_idx] = "X"                
        end
    end    
    return coords_data, coords_stab, string(types_stab...)
end

"""
    tensor_data(stab_type)

Return the tensor of certain type for the tensor network for the surface-hexagonal code
"""
function tensor_data(stab_type)
    @assert stab_type isa String

    stab_type == "Z" && return [0.0, 2.0]
    stab_type == "X" && return [0.0, 1.0]

    arr0 = tensor_data(stab_type[1:end-1])
    if stab_type[end] == 'X'
        arr1 = broadcast(mod, 1.0 .- arr0, 4.0)
    elseif stab_type[end] == 'Z'
        arr1 = broadcast(mod, arr0 .- 2.0, 4.0)
    end
    arr = cat(arr0, arr1, dims = length(size(arr0))+1)
    return arr
end

"""
    tn_template_surf_hex(d::Int)

Return the template for the tensor network for the surface-hexagonal GKP code at distance d. 

The template is setup in such a way that one only need to fill in all the entries of the tensors for the data qubits (with values 0,1,2,3), before contracting the tensor network

Assuming square surface code, so that we can use a single TN for all four probs [why???]

Returns:
    TN: The tempalte tensor network
    indices: The indices for the data qubits in the tensor network
"""
function tn_template_surf_hex(d::Int)
    dx, dz = d, d
    coords_data, coords_stab, types_stab = get_coords_surf_hex(dx, dz)
    Z_stabs = rectangular_surface_code_Z_stabilizers(dx, dz)
    X_stabs = rectangular_surface_code_X_stabilizers(dx, dz)
    stabilizers = Dict()
    for (_, stab) in sort(Z_stabs)
        stabilizers[length(stabilizers)+1] = stab
    end
    for (_, stab) in sort(X_stabs)
        stabilizers[length(stabilizers)+1] = stab
    end


    num_data_qubits = dx * dz
    num_stabilizers = Int(dx*dz-1)

    adj_data = Dict(k=>Int[] for k in 1 : num_data_qubits)
    for (key, val) in stabilizers
        for i in val
            push!(adj_data[i], key + num_data_qubits)
        end
    end

    TN = TensorNetwork(undef, 2*dx*dz-1);
    for i in 1 : num_data_qubits
        adj = adj_data[i]
        stab_type = types_stab[adj .- num_data_qubits]
#         println(i, stab_type, adj)        
        TN[i] = Tensor(adj, tensor_data(stab_type), coords_data[i][1], coords_data[i][2])
    end

    for i in 1 : num_stabilizers
        adj = stabilizers[i]

        if length(adj) == 2
            tensor = [x==y ? 1.0 : 0.0 for x=0:1, y=0:1]
        elseif length(adj) == 4
            tensor = [x==y&&y==z&&z==w ? 1.0 : 0.0 for x=0:1, y=0:1, z=0:1, w=0:1]
        end

        TN[i+num_data_qubits] = Tensor(adj, tensor, coords_stab[i][1], coords_stab[i][2])
    end

    SweepContractor.connect!(SweepContractor.hull!(TN))
    SweepContractor.sort!(TN)
    coords_TN = [(tn.x, tn.y) for tn in TN]

    indices = zeros(Int, num_data_qubits)
    for (ind1, coord) in enumerate(coords_data)
        ind2 = findfirst(x->norm((x[1], x[2]) .- coord) < 1e-10, coords_TN)
        indices[ind1] = ind2
    end

    return TN, indices
end


"""
    tn_surf_hex(ηs, σ, TN, indices, Z, X, χ; S = [2 1; 0 sqrt(3)] / (12)^(1/4), Nv=5)

Perform tensor network MLD for surface-hexagonal GKP code

Args:
    ηs: The candidate error in both the q and p subspaces
    σ: The noise strength
    TN: The template for the triangular color-square code
    indices: The indices for the data qubits in the tensor network
    Z: The Z logical operator of the triangular color code
    X: The X logical operator of the triangular color code
    χ: The max bond dimension
    S: The symplectic matrix for the inner one-mode GKP code
    Nv: The truncation for the number of closest points for the inner square GKP code

Returns:
    lstar: The most probable logical operator, either I or X
    prob_I: The coset probability for I
    prob_X: The coset probability for X
    prob_Y: The coset probability for Y
    prob_Z: The coset probability for Z
    num_truncs_I: The number of truncations in the tensor network when evaluating prob_I
    num_truncs_X: The number of truncations in the tensor network when evaluating prob_X
    num_truncs_Y: The number of truncations in the tensor network when evaluating prob_Y
    num_truncs_Z: The number of truncations in the tensor network when evaluating prob_Z
"""    
function tn_surf_hex(ηs, σ, TN, indices, Z, X, χ; S = [2 1; 0 sqrt(3)] / (12)^(1/4), Nv=5)
    TN0 = deepcopy(TN) # Backup of TN
    num_qubits = Int(length(ηs)/2)
    d = Int(sqrt(num_qubits))

    ws = get_ws_non_square(ηs, σ; S = S, Nv=Nv)

    TN = deepcopy(TN0)    
    for i in 1 : num_qubits
        for j in [0, 1, 2, 3] # I, X, Z, Y
            TN[indices[i]].arr[TN[indices[i]].arr .== j] .= ws[i][j+1]
        end
    end       
    
    # Get the intermediate result before evaluating the tensors for logicals
    MPS_t, MPS_i, resexp, num_truncs = sweep_contract_v2!(TN, χ, χ; 
                    start_ind=1, end_ind=length(TN)-d)
        
    # Get the prob for prob_I
    res_I, resexp_I, num_truncs_I = sweep_contract_v2!(TN, χ, χ; 
                    start_ind=length(TN)-d+1, MPS_t=deepcopy(MPS_t), MPS_i=deepcopy(MPS_i))
    
    num_truncs_I += num_truncs
    prob_I = log10(res_I) + (resexp_I + resexp) * log10(2)

    # Get the weights for the logical operators Z
    Z_logical = surface_code_Z_logicals(d)[1]
    TN[end-d+1:end] = deepcopy(TN0[end-d+1:end])

    # Modify the tensor network
    for i in Z_logical
        for j in [0, 1, 2, 3] # I, X, Z, Y => Z, Y, I, X
            TN[indices[i]].arr[TN[indices[i]].arr .== j] .= ws[i][mod(j+2, 4)+1]
        end
    end
    
    # Get the prob for prob_Z
    res_Z, resexp_Z, num_truncs_Z = sweep_contract_v2!(TN, χ, χ; 
                    start_ind=length(TN)-d+1, MPS_t=MPS_t, MPS_i=MPS_i)
    
    num_truncs_Z += num_truncs
    prob_Z = log10(res_Z) + (resexp_Z + resexp) * log10(2)
    

    ###### Get prob_X and prob_Y
    X_logical = surface_code_X_logicals(d)[1]
    TN = deepcopy(TN0)    
    for i in 1 : num_qubits
        for j in [0, 1, 2, 3] # I, X, Z, Y
            # If the i-th qubit is not part of X logical operator 
            # then j=0,1,2,3 => I, X, Z, Y
            # else j=0,1,2,3 => X, I, Y, Z
            i ∈ X_logical ? ws_val = ws[i][mod(5-j, 4)+1] : ws_val = ws[i][j+1]
            TN[indices[i]].arr[TN[indices[i]].arr .== j] .= ws_val
        end
    end

    # Get the intermediate result before evaluating the tensors for logicals
    MPS_t, MPS_i, resexp, num_truncs = sweep_contract_v2!(TN, χ, χ; 
                    start_ind=1, end_ind=length(TN)-d)
        
    # Get the prob for prob_X
    res_X, resexp_X, num_truncs_X = sweep_contract_v2!(TN, χ, χ; 
                    start_ind=length(TN)-d+1, MPS_t=deepcopy(MPS_t), MPS_i=deepcopy(MPS_i))
    
    num_truncs_X += num_truncs
    prob_X = log10(res_X) + (resexp_X + resexp) * log10(2)
    
    TN[end-d+1:end] = deepcopy(TN0[end-d+1:end])

    # Modify the tensor network
    for i in Z_logical
        for j in [0, 1, 2, 3] # I, X, Z, Y
            # For the first qubit, j=0,1,2,3 => Y, Z, X, I
            # else j=0,1,2,3 => Z, Y, I, X
            i == 1 ? ws_val = ws[i][3-j+1] : ws_val = ws[i][mod(j+2, 4)+1]
    
            TN[indices[i]].arr[TN[indices[i]].arr .== j] .= ws_val
        end
    end
    
    # Get the prob for prob_Z
    res_Y, resexp_Y, num_truncs_Y = sweep_contract_v2!(TN, χ, χ; 
                    start_ind=length(TN)-d+1, MPS_t=MPS_t, MPS_i=MPS_i)
    
    num_truncs_Y += num_truncs
    prob_Y = log10(res_Y) + (resexp_Y + resexp) * log10(2)

    #####    
    prob = max(prob_I, prob_X, prob_Y, prob_Z)
    if prob == prob_I
        lstar = zeros(length(ηs))
    elseif prob == prob_X
        lstar = X
    elseif prob == prob_Z
        lstar = Z
    elseif prob == prob_Y
        lstar = X+Z
    end
    return lstar, prob_I, prob_X, prob_Y, prob_Z, num_truncs_I, num_truncs_X, num_truncs_Y, num_truncs_Z
end


"""
    mwms_surface_square(ηs, σ, Kmax; subspace="z")

Decode surface-hexagonal GKP code using k-mwm method

Args:
    ηs: The candidate error in either the q or p subspace
    σ: The noise strength
    Kmax: The number of MWMs 
    subspace: The subspace in which decoding is performed, either "x" or "z"
    Nv: The number of closeset points for each dimension

Returns:
    p_I: Vector with Kmax+1 elements for the coset probability for I
    p_X: Vector with Kmax+1 elements for the coset probability for X

Notes: Here we calculate the accumulated coset probabilities as we increase k from 1 to Kmax.
    The (Kmax+1)-th coset probabilities are calculated by including all the candidate MWMs 
    found when identifying the Kmax MWMs.
"""    
function mwms_surface_square(ηs, σ, Kmax; subspace="z", Nv=5)

    d = √(length(ηs))
    if !isinteger(d)
        error("The length of ηs has to be d^2 for an integer d")
    end
    d = Int(d)
    
    @assert subspace ∈ ["x", "z"]
    
    if subspace=="z"
        stabilizers = surface_code_Z_stabilizers(d)
        logicals = surface_code_Z_logicals(d)
    else
        stabilizers = surface_code_X_stabilizers(d)
        logicals = surface_code_X_logicals(d)
    end
    
    g, hvs, edge_to_mode_dict, closest_integers, second_closest_integers, G, vertices_edge_mapping = decoding_graph(ηs, stabilizers)
    error_list = mod.(closest_integers, 2)

    Matchings, Matchings2 = mwms(g, hvs, G, vertices_edge_mapping, Kmax)

    if hvs == zeros(length(hvs))
        pushfirst!(Matchings, [])
        Matchings = Matchings[1:end-1]
    end
    
    function get_p_I_p_X_from_path(p_I, p_X, path)
        correction_list = zeros(Int, length(ηs))
        for (i, j) in path
            if Set([i,j]) ∈ keys(edge_to_mode_dict)
                correction_list[edge_to_mode_dict[Set([i,j])][1]] += 1
            end
        end

        cp = [
            if (correction == 1) second_closest_integers[i] else closest_integers[i] end 
                for (i, correction) in enumerate(correction_list)
        ] .* √π
        cps = (2*√π) .* closest_points_Zn((ηs - cp)./(2*√π), Nv)
        cps = [cp2 + cp for cp2 in cps]
        prob = sum([gaussian(σ, ηs - cp2) for cp2 in cps])

        final_list = correction_list - error_list
        nx = mod(sum(final_list[logicals[1]]), 2)
        mod(nx, 2) == 0 ? (p_I += prob) : (p_X += prob)
        return p_I, p_X
    end

    p_I, p_X = 0, 0
    ps_I = []
    ps_X = []
    for (ind, path) in enumerate(Matchings)
        p_I, p_X = get_p_I_p_X_from_path(p_I, p_X, path)
        push!(ps_I, p_I)
        push!(ps_X, p_X)
    end

    for (ind, path) in enumerate(Matchings2)
        p_I, p_X = get_p_I_p_X_from_path(p_I, p_X, path)
    end
    push!(ps_I, p_I)
    push!(ps_X, p_X)
    return ps_I, ps_X
end


"""
    mwms_surface_square(ηs, σ, Kmax; subspace="z")

Decode surface-hexagonal GKP code using k-mwm method

Args:
    ηs: The candidate error in either the q or p subspace
    σ: The noise strength
    Kmax: The number of MWMs 
    subspace: The subspace in which decoding is performed, either "x" or "z"
    Nv: The number of closeset points for each dimension

Returns:
    p_I: Vector with Kmax+1 elements for the coset probability for I
    p_X: Vector with Kmax+1 elements for the coset probability for X

Notes: Here we calculate the accumulated coset probabilities as we increase k from 1 to Kmax.
    The (Kmax+1)-th coset probabilities are calculated by including all the candidate MWMs 
    found when identifying the Kmax MWMs.
"""    
function mwms_surface_square_v2(ηs, σ, Kmax; subspace="z", Nv=5)

    d = √(length(ηs))
    if !isinteger(d)
        error("The length of ηs has to be d^2 for an integer d")
    end
    d = Int(d)
    
    @assert subspace ∈ ["x", "z"]
    
    if subspace=="z"
        stabilizers = surface_code_Z_stabilizers(d)
        logicals = surface_code_Z_logicals(d)
    else
        stabilizers = surface_code_X_stabilizers(d)
        logicals = surface_code_X_logicals(d)
    end
    
    # g, hvs, edge_to_mode_dict, closest_integers, second_closest_integers, G, vertices_edge_mapping = decoding_graph(ηs, stabilizers)
    g, hvs, edge_to_mode_dict, closest_integers, second_closest_integers, G, vertices_edge_mapping, edge_vertices_mapping = decoding_graph_v2(ηs, stabilizers)

    error_list = mod.(closest_integers, 2)

    # Matchings, Matchings2 = mwms(g, hvs, G, vertices_edge_mapping, Kmax)
    Matchings, Matchings2 = mwms_v2(g, hvs, G, vertices_edge_mapping, edge_vertices_mapping, Kmax)

    if hvs == zeros(length(hvs))
        pushfirst!(Matchings, [])
        Matchings = Matchings[1:end-1]
    end
    
    function get_p_I_p_X_from_path(p_I, p_X, path)
        correction_list = zeros(Int, length(ηs))
        for e in path
            (i, j) = vertices_edge_mapping[e]
            if Set([i,j]) ∈ keys(edge_to_mode_dict)
                correction_list[edge_to_mode_dict[Set([i,j])][1]] += 1
            end
        end

        cp = [
            if (correction == 1) second_closest_integers[i] else closest_integers[i] end 
                for (i, correction) in enumerate(correction_list)
        ] .* √π
        cps = (2*√π) .* closest_points_Zn((ηs - cp)./(2*√π), Nv)
        cps = [cp2 + cp for cp2 in cps]
        prob = sum([gaussian(σ, ηs - cp2) for cp2 in cps])

        final_list = correction_list - error_list
        nx = mod(sum(final_list[logicals[1]]), 2)
        mod(nx, 2) == 0 ? (p_I += prob) : (p_X += prob)
        return p_I, p_X
    end

    p_I, p_X = 0, 0
    ps_I = []
    ps_X = []
    for (ind, path) in enumerate(Matchings)
        p_I, p_X = get_p_I_p_X_from_path(p_I, p_X, path)
        push!(ps_I, p_I)
        push!(ps_X, p_X)
    end

    for (ind, path) in enumerate(Matchings2)
        p_I, p_X = get_p_I_p_X_from_path(p_I, p_X, path)
    end
    push!(ps_I, p_I)
    push!(ps_X, p_X)
    return ps_I, ps_X
end


"""
    mwms_surface_square(ηs, σ, Kmax; S_hex = [2 1; 0 √3] / (12)^(1/4), Nv=5)

Decode surface-hexagonal GKP code using k-mwm method

Args:
    ηs: The candidate error in either the q or p subspace
    σ: The noise strength
    Kmax: The number of MWMs 
    S_hex: The symplectic matrix for the one mode hexagonal code
    Nv: The number of closeset points for each dimension

Returns:
    ps_I: Vector with Kmax elements for the coset probability for I
    ps_X: Vector with Kmax elements for the coset probability for X
    ps_Y: Vector with Kmax elements for the coset probability for Y
    ps_Z: Vector with Kmax elements for the coset probability for Z
"""
function mwms_surface_hexagonal(ηs, σ, Kmax; S_hex = [2 1; 0 √3] / (12)^(1/4), Nv=5)
    d = √(length(ηs)/2)
    if !isinteger(d)
        error("The length of ηs has to be 2d^2 for an integer d")
    end
    d = Int(d)    
    surface_code_x_stabilizers = surface_code_X_stabilizers(d)
    surface_code_z_stabilizers = surface_code_Z_stabilizers(d)

    S_hex_T = Matrix(transpose(S_hex))

    M0 = surface_code_M(d)
    M = M0 * BlockDiagonal([S_hex_T for _ in 1 : d^2])

    Mperp = GKP_logical_operator_generator(M)
    Ω = Ω_matrix(Mperp)
    transposeMperp_qq = transpose(Mperp)[1:2:end, 1:2:end]
    transposeMperp_pp = transpose(Mperp)[2:2:end, 2:2:end]    
    
    function get_cp_from_a_matching(
            matching, 
            num_mode, 
            vertices_edge_mapping, 
            edge_to_mode_dict,
            closest_integers,
            second_closest_integers,
            Γ
        )
        correction_list = zeros(Int, num_mode)

        for e in matching
            (i, j) = vertices_edge_mapping[e]
            if Set([i,j]) ∈ keys(edge_to_mode_dict)
                correction_list[edge_to_mode_dict[Set([i,j])][1]] += 1
            end
        end

        cp = [
            if (correction == 1) second_closest_integers[i] else closest_integers[i] end 
                for (i, correction) in enumerate(correction_list)
        ] .* √π .* Γ
        return cp
    end
    
    function decode_subspace(ηs2, surface_code_stabilizers, Γ, Kmax, transposeMperp, subspace)
        g, hvs, edge_to_mode_dict, closest_integers, second_closest_integers, G, vertices_edge_mapping, edge_vertices_mapping = decoding_graph_v2(ηs2, surface_code_stabilizers, Γ)

        if subspace == "p"
            Matchings, _ = mwms_v2(g, hvs, G, vertices_edge_mapping, edge_vertices_mapping, Kmax)
        elseif subspace == "q"
            Matchings, X_q = mwms_v2(g, hvs, G, vertices_edge_mapping, edge_vertices_mapping, Kmax; option=2)
        end
        
        cps2 = [
            get_cp_from_a_matching(
                matching,                 
                length(ηs2), 
                vertices_edge_mapping, 
                edge_to_mode_dict,
                closest_integers,
                second_closest_integers,
                Γ
            )
            for matching in Matchings
        ]

        as_2 = []
        for cp in cps2
            u = inv(√(2π) * transposeMperp) * cp
            @assert round.(Int, u) ≈ u
            push!(as_2, round.(Int, u))
        end

        if subspace == "p"
            return as_2
        elseif subspace == "q"
            return as_2, g, edge_to_mode_dict, closest_integers, second_closest_integers, G, vertices_edge_mapping, edge_vertices_mapping, Matchings, X_q
        end
    end

    ηsp = ηs[2:2:end]
    Γ4 = S_hex[2,2] * ones(length(ηsp))

    ηsq = ηs[1:2:end]
    Γ1 = S_hex[1,1] * ones(length(ηsq))
    Γ2 = S_hex[1,2] * ones(length(ηsq))        
    
    as_p = decode_subspace(ηsp, surface_code_x_stabilizers, Γ4, Kmax, transposeMperp_pp, "p")

    g_q_list = []
    edge_to_mode_dict_q_list = []
    closest_integers_q_list = []
    second_closest_integers_q_list = []
    G_q_list = []
    vertices_edge_mapping_q_list = []
    edge_vertices_mapping_q_list = []
    XK_q_list = []
    X_q_list = []
    X = PriorityQueue{Int, Real}()
    for (ind_p, a_p) in enumerate(as_p)
        ηsq2 = ηsq - √(2π) * Γ2 .* (inv(M0[1:2:end, 1:2:end]) * a_p)

        a_q, g_q, edge_to_mode_dict_q, closest_integers_q, second_closest_integers_q, G_q, vertices_edge_mapping_q, edge_vertices_mapping_q, XK_q, X_q = decode_subspace(ηsq2, surface_code_z_stabilizers, Γ1, 1, transposeMperp_qq, "q")
        a_q = a_q[1]

        push!(g_q_list, g_q)
        push!(edge_to_mode_dict_q_list, edge_to_mode_dict_q)
        push!(closest_integers_q_list, closest_integers_q)
        push!(second_closest_integers_q_list, second_closest_integers_q)    
        push!(G_q_list, G_q)
        push!(vertices_edge_mapping_q_list, vertices_edge_mapping_q)
        push!(edge_vertices_mapping_q_list, edge_vertices_mapping_q)
        push!(XK_q_list, XK_q)    
        push!(X_q_list, X_q)

        a = zeros(Int, 2length(a_q))
        a[1:2:end] .= round.(Int, a_q)
        a[2:2:end] .= round.(Int, a_p)

        y = √(2π) * transpose(Mperp) * a
        X[ind_p] = norm(y - ηs)
    end    

    as = []
    for k in 1 : Kmax
        ind_q = dequeue!(X)
        g = g_q_list[ind_q]
        edge_to_mode_dict = edge_to_mode_dict_q_list[ind_q]
        closest_integers = closest_integers_q_list[ind_q]
        second_closest_integers = second_closest_integers_q_list[ind_q]    
        G = G_q_list[ind_q]
        vertices_edge_mapping = vertices_edge_mapping_q_list[ind_q]
        edge_vertices_mapping = edge_vertices_mapping_q_list[ind_q]
        XK = XK_q_list[ind_q]
        X_q = X_q_list[ind_q]
            
        Mnext = XK[end]
        cp = get_cp_from_a_matching(
                Mnext,                 
                length(closest_integers), 
                vertices_edge_mapping, 
                edge_to_mode_dict,
                closest_integers,
                second_closest_integers,
                Γ1
            )

        u = inv(√(2π) * transposeMperp_qq) * cp
        @assert round.(Int, u) ≈ u

        a = zeros(Int, 2length(u))
        a[1:2:end] .= round.(Int, u)
        a[2:2:end] .= round.(Int, as_p[ind_q])    

        push!(as, a)
            
        X_q = get_next_matching!(g, X_q, XK, G, vertices_edge_mapping, edge_vertices_mapping)

        if length(X_q) == 0
            # Discard X[ind_q] if it ends
            continue
        end
        
        (Mnext, _, _, _), w_Mnext = peek(X_q)
        if w_Mnext == Inf
            # Discard X[ind_q]
            continue
        end
            
        y = √(2π) * transpose(Mperp) * a
        X[ind_q] = norm(y - ηs) 

        push!(XK_q_list[ind_q], Mnext)
        X_q_list[ind_q] = X_q
    end

    # Get p_I, p_X, p_Y, p_Z
    X_logical = surface_code_X_logicals(d)[1]
    X = zeros(2d^2)
    X[2 .* X_logical .- 1] .= 1
    Z_logical = surface_code_Z_logicals(d)[1]
    Z = zeros(2d^2)
    Z[2 .* Z_logical] .= 1

    Z = BlockDiagonal([S_hex for _ in 1:d^2]) * Z .* √π
    X = BlockDiagonal([S_hex for _ in 1:d^2]) * X .* √π

    p_I, p_X, p_Y, p_Z = 0, 0, 0, 0
    ps_I, ps_X, ps_Y, ps_Z = [], [], [], []
    for a in as
        y = √(2π) * transpose(Mperp) * a
        res_error = y-ηs
        
        # Get the probability of res_error
        prob = 1
        for i in 1 : Int(length(res_error)/2)
            χi = res_error[2i-1:2i]/√π
            cps, _ = closest_points(χi, 2 * S_hex_T, Nv)
            prob *= sum([exp(-π * norm(cp-χi)^2) for cp in cps])
        end        
        
        nx = abs(transpose(y) * Ω * Z / π)
        nz = abs(transpose(y) * Ω * X / π)
        nx = mod(round(nx), 2)
        nz = mod(round(nz), 2)    

        @assert abs(round(Int, nx)-nx)<1e-5
        @assert abs(round(Int, nz)-nz)<1e-5

        if nx ≈ 0 && nz ≈ 0
            p_I += prob
        elseif nx ≈ 0 && nz ≈ 1
            p_Z += prob
        elseif nx ≈ 1 && nz ≈ 0
            p_X += prob
        elseif nx ≈ 1 && nz ≈ 1
            p_Y += prob    
        end
        push!(ps_I, p_I)
        push!(ps_X, p_X)
        push!(ps_Y, p_Y)
        push!(ps_Z, p_Z)
    end
    
    return ps_I, ps_X, ps_Y, ps_Z
end