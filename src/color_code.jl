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

## Triangular color codes

"""
    triangular_color_code_num_qubits(d::Int)

Return the numbef of modes in the triangular color code.
"""
function triangular_color_code_num_qubits(d::Int)
    if mod(d, 2) != 1 || d < 1
        error("The distance of a color code should be an odd number bigger or equal to 1.")
    else
        return Int((3d^2+1)/4)    
    end 
end

"""
    triangular_color_code_stabilizers(d::Int, color::String)

Return the stabilizers for the triangular color code at distance d with certain color ("R", "G", "B")
"""
function triangular_color_code_stabilizers(d::Int, color::String)

    if !(color ∈ ("R", "G", "B"))
        error("`color` has to be either `R`, `G` or `B`.")
    end    

    if d==3
        dict = Dict{Int64, Vector{Int64}}()
        if color=="G"
            dict[1] = [1,2,3,4] # Green
        elseif color=="B"
            dict[1] = [3,4,6,7] # Blue
        else
            dict[1] = [2,3,5,6] # Red
        end
        return dict
    end

    N_old = triangular_color_code_num_qubits(d-2)

    dict = triangular_color_code_stabilizers(d-2, color)
    num_old_stabs = length(dict)

    if color=="B"
        # We first complete the last (d-3)/2 blue stabilizers from weight 4 to weight 6
        for i in 1 : Int((d-3)/2)
            old_stab = dict[num_old_stabs-(i-1)]
            dict[num_old_stabs-(i-1)] = [old_stab..., old_stab[end]+d-3, old_stab[end]+d-2]
        end    
    end

    for i in 1 : Int((d-1)/2)

        if color=="G"

            # Add the green stabilizer
            if i == Int((d-1)/2) # Boundary weight 4 green stabilizer
                dict[num_old_stabs+i] = [
                    N_old, N_old + d-2, N_old + 2(d-2), N_old + 2(d-2)+1
                ]
            else
                dict[num_old_stabs+i] = [
                    N_old-(d-2) + 2i-1, N_old-(d-2) + 2i-1+1, 
                    N_old + 2i-1, N_old + 2i, 
                    N_old + 2i-1+(d-2), N_old + 2i+(d-2), 
                ]
            end
        elseif color=="R"

            # Add the red stabilizer
            if i == 1 # Boundary weight 4 red stabilizer
                dict[num_old_stabs+i] = [
                    N_old+1, N_old+1 + d-2, N_old+1 + d-2 + d-1, N_old+1 + d-2 + d
                ]
            else
                dict[num_old_stabs+i] = [
                    N_old+2 + 2i-4, N_old+2 + 2i-4+1, 
                    N_old+2 + 2i-4+(d-2), N_old+2 + 2i-4+1+(d-2), 
                    N_old+2 + 2i-4+(d-2)+d, N_old+2 + 2i-4+(d-2)+d+1
                ]
            end
        else    
            # Add the blue stabilizer
            dict[num_old_stabs+i] = [
                N_old+(d-2) + 2i-1, N_old+(d-2) + 2i, N_old+(d-2) + 2i-1+d, N_old+(d-2) + 2i+d
            ]
        end
    end

    return dict
end

"""
    triangular_color_code_stabilizers(d::Int)

Return the stabilizers for the triangular color code at distance d
"""
function triangular_color_code_stabilizers(d::Int)
    if mod(d, 2) != 1 || d < 3
        error("The distance of a color code should be an odd number bigger or equal to 3.")
    end
    
    if d==3
        X_dict = Dict{Int64, Vector{Int64}}()
        X_dict[1] = [1,2,3,4] # Green
        X_dict[2] = [2,3,5,6] # Red
        X_dict[3] = [3,4,6,7] # Blue
        return X_dict
    end
    
    X_dict = triangular_color_code_stabilizers(d-2)
    num_old_stabs = length(X_dict)
    
    # We first complete the last (d-3)/2 blue stabilizers from weight 4 to weight 6
    for i in 1 : Int((d-3)/2)
        old_stab = X_dict[num_old_stabs-3(i-1)]
        X_dict[num_old_stabs-3(i-1)] = [old_stab..., old_stab[end]+d-3, old_stab[end]+d-2]
    end
    
    N_old = Int((3(d-2)^2+1)/4)
    for i in 1 : Int((d-1)/2)
        # Add the green stabilizer
        if i == Int((d-1)/2) # Boundary weight 4 green stabilizer
            X_dict[num_old_stabs+3i-2] = [
                N_old, N_old + d-2, N_old + 2(d-2), N_old + 2(d-2)+1
            ]
        else
            X_dict[num_old_stabs+3i-2] = [
                N_old-(d-2) + 2i-1, N_old-(d-2) + 2i-1+1, 
                N_old + 2i-1, N_old + 2i, 
                N_old + 2i-1+(d-2), N_old + 2i+(d-2), 
            ]
        end
        
        # Add the red stabilizer
        if i == 1 # Boundary weight 4 red stabilizer
            X_dict[num_old_stabs+3i-1] = [
                N_old+1, N_old+1 + d-2, N_old+1 + d-2 + d-1, N_old+1 + d-2 + d
            ]
        else
            X_dict[num_old_stabs+3i-1] = [
                N_old+2 + 2i-4, N_old+2 + 2i-4+1, 
                N_old+2 + 2i-4+(d-2), N_old+2 + 2i-4+1+(d-2), 
                N_old+2 + 2i-4+(d-2)+d, N_old+2 + 2i-4+(d-2)+d+1
            ]
        end
                
        # Add the blue stabilizer
        X_dict[num_old_stabs+3i] = [
            N_old+(d-2) + 2i-1, N_old+(d-2) + 2i, N_old+(d-2) + 2i-1+d, N_old+(d-2) + 2i+d
        ]
    end
    
    return X_dict
end

"""
    merge_two_stabilizers(stabilizers_1::Dict{Int64, Vector{Int64}}, stabilizers_2::Dict{Int64, Vector{Int64}})

Merge two stabilizers of type Dict{Int64, Vector{Int64}}
"""
function merge_two_stabilizers(
    stabilizers_1::Dict{Int64, Vector{Int64}}, 
    stabilizers_2::Dict{Int64, Vector{Int64}}
)
    stabilizers = Dict{Int64, Vector{Int64}}()
    for (_, val) in sort(stabilizers_1)
        stabilizers[length(stabilizers)+1] = val
    end

    for (_, val) in sort(stabilizers_2)
        stabilizers[length(stabilizers)+1] = val
    end
    return stabilizers
end


"""
    triangular_color_code_logicals(d::Int)

Return the logical operator for the triangular color code at distance d
"""
function triangular_color_code_logicals(d::Int)
    if d==3
        dict = Dict{Int64, Vector{Int64}}()
        dict[1] = [1, 2, 5]
        return dict
    end    
    N_old = triangular_color_code_num_qubits(d-2)
    N_new = triangular_color_code_num_qubits(d)    
    dict = triangular_color_code_logicals(d-2)
    dict[1] = [dict[1]..., N_old+1, N_new-(d-1)]

    return dict
end

"""
    triangular_color_code_normalizers(d::Int)

Return the normalizers (stabilizers and logical operator) for the triangular color code at distance d
"""
function triangular_color_code_normalizers(d::Int)
    dict = triangular_color_code_stabilizers(d)
    dict[length(dict)+1] = triangular_color_code_logicals(d)[1]
    return dict
end

"""
    triangular_color_code_Mq(d::Int)

Return the generator for the color code of distance d in the q subspace
"""
function triangular_color_code_Mq(d::Int)
    Mq = 2 * diagm(ones(Int64, triangular_color_code_num_qubits(d)))
    dict = triangular_color_code_stabilizers(d)
    for (_, value_list) in dict
        min_value = minimum(value_list)
        for value in value_list
            Mq[min_value, value] = 1
        end
    end
    Mq = Mq/√2 # To make sure det(Mq) = √2
    return Mq
end

"""
    triangular_color_code_Mp(d::Int)

Return the generator for the triangular color code of distance d in the p subspace
"""
triangular_color_code_Mp(d::Int) = triangular_color_code_Mq(d)

"""
    triangular_color_code_M(d::Int)

Return the generator for the triangular color code of distance d
"""
function triangular_color_code_M(d::Int)
    Mq = triangular_color_code_Mq(d)
    M = BlockDiagonal([Mq, Mq])
    T = basis_transformation(size(Mq)[1])
    return transpose(T) * M * T
end

######## functions for the tensor network MLD for the color code

"""
    get_coords_triangular_color_codes(d::Int; shift=0.001)

Return the coordinates for the tensor network for the triangular color code at distance d. The shift is for shifting the coordinates for the tensors of the stabilizers a bit so that the tensor network can be contracted properly when using the SweepContractor.jl package.
"""
function get_coords_triangular_color_codes(d::Int; shift=0.001)
    if mod(d, 2) != 1 || d < 3
        error("The distance of a color code should be an odd number bigger or equal to 3.")
    end

    R = [cos(pi/3) sin(pi/3); -sin(pi/3) cos(pi/3)]

    if d==3
        coords_data = [(0, 0), (-1/2, -√3/2), (0, -√3), (1, -√3), (-3/2, -3√3/2), (-1/2, -3√3/2), (3/2, -3√3/2)]
        coords_stab = [(1/2, -√3/2), (-1, -√3), (1/2, -3√3/2)]
    else
        coords_data, coords_stab = get_coords_triangular_color_codes(d-2)
        coords_data = Tuple.([inv(R) * collect(coord) for coord in coords_data])
        coords_stab = Tuple.([inv(R) * collect(coord) for coord in coords_stab])

        # Add three layers of data qubits
        coords_data_1 = coords_data[end-(d-2)+1:end]
        coords_data_2 = coords_data[end-2(d-2)+2:end-(d-2)]

        pushfirst!(coords_data_2, (coords_data_2[1][1]-2, coords_data_2[1][2]))
        coords_data_2 = [(coord[1], coord[2]-√3) for coord in coords_data_2]

        coords_data_1 = [(coord[1], coord[2]-√3) for coord in coords_data_1]
        push!(coords_data_1, (coords_data_1[end][1]+1, coords_data_1[end][2]))

        coords_data_3 = [(coord[1], coord[2]-√3) for coord in coords_data_2]
        pushfirst!(coords_data_3, (coords_data_3[1][1]-1, coords_data_3[1][2]))
        push!(coords_data_3, (coords_data_3[end][1]+2, coords_data_3[1][2]))

        coords_data = vcat(coords_data, coords_data_2, coords_data_1, coords_data_3)

        # Add the stabilizers
        coords_stab_old = coords_stab[end - 3 * Int((d-3)/2)+1 : end]
        coords_stab_old = [(coord[1] - 3/2, coord[2]-3√3/2) for coord in coords_stab_old]
        coords_stab_old = vcat(coords_stab_old, [(coord[1] +3, coord[2]) for coord in coords_stab_old[end-2:end]]...)
        coords_stab = vcat(coords_stab, coords_stab_old...)        
    end

    coords_data = Tuple.([R * collect(coord) for coord in coords_data])
    coords_stab = Tuple.([R * collect(coord) for coord in coords_stab])
    coords_stab = [coord .- (0, shift) for coord in coords_stab]
    return coords_data, coords_stab
end

"""
    tn_template_color_square(d::Int)

Return the template for the tensor network for the color-square GKP code at distance d. 

The template is setup in such a way that one only need to fill in all the entries of the tensors for the data qubits (with values -1 and 1), before contracting the tensor network

Returns:
    TN: The tempalte tensor network
    indices: The indices for the data qubits in the tensor network
"""
function tn_template_color_square(d::Int)
    coords_data, coords_stab = get_coords_triangular_color_codes(d)
    stabilizers = triangular_color_code_stabilizers(d)

    num_data_qubits = Int((3d^2+1)/4)
    num_stabilizers = Int((num_data_qubits-1)/2)

    adj_data = Dict(k=>Int[] for k in 1 : num_data_qubits)
    for (key, val) in stabilizers
        for i in val
            push!(adj_data[i], key + num_data_qubits)
        end
    end

    TN = TensorNetwork(undef, num_data_qubits + num_stabilizers);

    for i in 1 : num_data_qubits
        adj = adj_data[i]

        # The "-1" slots in the tensors for the data qubits are to be 
        # filled before the tensor network is contracted.
        if length(adj) == 1
            tensor = [mod(x, 2)==1 ? -1.0 : 1.0 for x=0:1]
        elseif length(adj) == 2
            tensor = [mod(x+y, 2)==1 ? -1.0 : 1.0 for x=0:1, y=0:1]
        elseif length(adj) == 3
            tensor = [mod(x+y+z, 2)==1 ? -1.0 : 1.0 for x=0:1, y=0:1, z=0:1]
        end
        TN[i] = Tensor(adj, tensor, coords_data[i][1], coords_data[i][2])
    end

    for i in 1 : num_stabilizers
        adj = stabilizers[i]
        if length(adj) == 4
            tensor = [x==y&&y==z&&z==w ? 1.0 : 0.0 for x=0:1, y=0:1, z=0:1, w=0:1]
        elseif length(adj) == 6
            tensor = [x1==x2&&x2==x3&&x3==x4&&x4==x5&&x5==x6 ? 1.0 : 0.0 for x1=0:1, x2=0:1, x3=0:1, x4=0:1, x5=0:1, x6=0:1]
        end    
        TN[i+num_data_qubits] = Tensor(adj, tensor, coords_stab[i][1], coords_stab[i][2])
    end
    
    # coords_data_2 = [(coord[1], coord[2]) for coord in coords_data]
    indices = zeros(Int, num_data_qubits)
    
    SweepContractor.connect!(SweepContractor.hull!(TN))
    SweepContractor.sort!(TN)

    coords_TN = [(tn.x, tn.y) for tn in TN]

    for (ind1, coord) in enumerate(coords_data)
        ind2 = findfirst(x->norm((x[1], x[2]) .- coord) < 1e-10, coords_TN)
        indices[ind1] = ind2
    end
        
    return TN, indices
end

"""
    tn_color_square(ηs_q, σ, logical, TN, indices, χ; Nv=5)

Perform tensor network MLD for color-square GKP code

Args:
    ηs_q: The candidate error in the q subspace
    σ: The noise strength
    logical: The logical operator of the triangular color code
    TN: The template for the triangular color-square code
    indices: The indices for the data qubits in the tensor network
    χ: The max bond dimension
    Nv: The truncation for the number of closest points for the inner square GKP code

Returns:
    lstar: The most probable logical operator, either I or X
    prob_I: The coset probability for I
    prob_X: The coset probability for X
    num_truncs_I: The number of truncations in the tensor network when evaluating prob_I
    num_truncs_X: The number of truncations in the tensor network when evaluating prob_X
"""
function tn_color_square(ηs_q, σ, logical, TN, indices, χ; Nv=5)
    num_qubits = length(ηs_q)
    d = Int(sqrt((4 * num_qubits - 1)/3)) # (3d^2+1)/4
    X = zeros(num_qubits)
    X[logical] .= 1/√2 * √(2π)
    TN0 = deepcopy(TN)
    
    TN = deepcopy(TN0)

    ws = get_ws_square(ηs_q, σ, Nv=Nv)
    w2s = get_ws_square(ηs_q .+ √π, σ, Nv=Nv)
    
    for i in 1 : num_qubits
        TN[indices[i]].arr[TN[indices[i]].arr .== -1] .= ws[i]
        TN[indices[i]].arr[TN[indices[i]].arr .== 1] .= w2s[i]
    end
            
    # Get the intermediate result before evaluating the tensors for logicals
    MPS_t, MPS_i, resexp, num_truncs = sweep_contract_v2!(TN, χ, χ; 
                    start_ind=1, end_ind=length(TN)-d)
        
    # Get the prob for prob_I
    res_I, resexp_I, num_truncs_I = sweep_contract_v2!(TN, χ, χ; 
                    start_ind=length(TN)-d+1, MPS_t=deepcopy(MPS_t), MPS_i=deepcopy(MPS_i))
    
    num_truncs_I += num_truncs
    prob_I = log10(res_I) + (resexp_I + resexp) * log10(2)
    
    # Modify the tensor network
    TN[end-d+1:end] = deepcopy(TN0[end-d+1:end])
    for i in logical
        TN[indices[i]].arr[TN[indices[i]].arr .== -1] .= w2s[i]
        TN[indices[i]].arr[TN[indices[i]].arr .== 1] .= ws[i]
    end        
    
    # Get the prob for prob_X
    res_X, resexp_X, num_truncs_X = sweep_contract_v2!(TN, χ, χ; 
                    start_ind=length(TN)-d+1, MPS_t=deepcopy(MPS_t), MPS_i=deepcopy(MPS_i))
    
    num_truncs_X += num_truncs
    prob_X = log10(res_X) + (resexp_X + resexp) * log10(2)
        
    prob_I < prob_X ? lstar = X : lstar = zeros(length(X))
    return lstar, prob_I, prob_X, num_truncs_I, num_truncs_X
end

"""
    tn_template_color_hex(d::Int)

Return the template for the tensor network for the color-hexagonal GKP code at distance d. 

The template is setup in such a way that one only need to fill in all the entries of the tensors for the data qubits (with values 0,1,2,3), before contracting the tensor network

Returns:
    TN: The tempalte tensor network
    indices: The indices for the data qubits in the tensor network    
"""
function tn_template_color_hex(d::Int; shift::Float64 = 0.001)
    coords_data, coords_stab = get_coords_triangular_color_codes(d; shift=shift)

    stabilizers = triangular_color_code_stabilizers(d)

    num_data_qubits = Int((3d^2+1)/4)
    num_stabilizers_z = Int((num_data_qubits-1)/2)

    adj_data = Dict(k=>Int[] for k in 1 : num_data_qubits)
    for (key, val) in stabilizers
        for i in val
            push!(adj_data[i], key + num_data_qubits)
        end
    end

    TN = LabelledTensorNetwork{Int}() ;

    for i in 1 : num_data_qubits
        adj = adj_data[i]

        if length(adj) == 1
            tensor = [x for x=0:3]
        elseif length(adj) == 2
            tensor = [xor(x,y) for x=0:3, y=0:3]
        elseif length(adj) == 3
            tensor = [xor(x,xor(y,z)) for x=0:3, y=0:3, z=0:3]
        end        
        TN[i] = Tensor(adj, tensor, coords_data[i][1], coords_data[i][2])
    end

    for i in 1 : num_stabilizers_z
        adj = stabilizers[i]
        if length(stabilizers[i]) == 4
            tensor = [x==y&&y==z&&z==w ? 1.0 : 0.0 for x=0:3, y=0:3, z=0:3, w=0:3]
        elseif length(stabilizers[i]) == 6
            tensor = [x1==x2&&x2==x3&&x3==x4&&x4==x5&&x5==x6 ? 1.0 : 0.0 for x1=0:3, x2=0:3, x3=0:3, x4=0:3, x5=0:3, x6=0:3]
        end    

        TN[i+num_data_qubits] = Tensor(adj, tensor, coords_stab[i][1], coords_stab[i][2])
    end

    TN = SweepContractor.delabel(TN)
    SweepContractor.planarise!(TN)
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
    tn_color_hex(ηs, σ, TN, indices, Z, X, χ; S = [2 1; 0 sqrt(3)] / (12)^(1/4), Nv=5)

Perform tensor network MLD for color-hexagonal GKP code

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
function tn_color_hex(ηs, σ, TN, indices, Z, X, χ; S = [2 1; 0 sqrt(3)] / (12)^(1/4), Nv=5)
    TN0 = deepcopy(TN) # Backup of TN
    num_qubits = Int(length(ηs)/2)
    d = Int(sqrt((4 * num_qubits - 1)/3)) # (3d^2+1)/4

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
    logical = triangular_color_code_logicals(d)[1]
    TN[end-d+1:end] = deepcopy(TN0[end-d+1:end])

    # Modify the tensor network
    for i in logical
        for j in [0, 1, 2, 3] # I, X, Z, Y => Z, Y, I, X
            TN[indices[i]].arr[TN[indices[i]].arr .== j] .= ws[i][mod(j+2, 4)+1]
        end
    end
    
    # Get the prob for prob_Z
    res_Z, resexp_Z, num_truncs_Z = sweep_contract_v2!(TN, χ, χ; 
                    start_ind=length(TN)-d+1, MPS_t=deepcopy(MPS_t), MPS_i=deepcopy(MPS_i))
    
    num_truncs_Z += num_truncs
    prob_Z = log10(res_Z) + (resexp_Z + resexp) * log10(2)
    
    # Get the weights for the logical operators X
    TN[end-d+1:end] = deepcopy(TN0[end-d+1:end])

    # Modify the tensor network
    for i in logical
        for j in [0, 1, 2, 3] # I, X, Z, Y => X, I, Y, Z
            TN[indices[i]].arr[TN[indices[i]].arr .== j] .= ws[i][mod(5-j, 4)+1]
        end
    end
    
    # Get the prob for prob_X
    res_X, resexp_X, num_truncs_X = sweep_contract_v2!(TN, χ, χ; 
                    start_ind=length(TN)-d+1, MPS_t=deepcopy(MPS_t), MPS_i=deepcopy(MPS_i))
    
    num_truncs_X += num_truncs
    prob_X = log10(res_X) + (resexp_X + resexp) * log10(2)
    
    # Get the weights for the logical operators Y
    TN[end-d+1:end] = deepcopy(TN0[end-d+1:end])

    # Modify the tensor network
    for i in logical
        for j in [0, 1, 2, 3] # I, X, Z, Y => Y, Z, X, I
            TN[indices[i]].arr[TN[indices[i]].arr .== j] .= ws[i][3-j+1]
        end
    end
    
    # Get the prob for prob_Y
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
#         println([num_truncs_I, num_truncs_X, num_truncs_Y, num_truncs_Z])
    return lstar, prob_I, prob_X, prob_Y, prob_Z, num_truncs_I, num_truncs_X, num_truncs_Y, num_truncs_Z
end



#################### Octagonal color codes

"""
    octagonal_color_code_num_qubits(d::Int)

Return the numbef of modes in the octagonal color code.
"""
function octagonal_color_code_num_qubits(d::Int)
    if mod(d, 2) != 0 || d < 2
        error("The distance of a periodic octagon color code should be an even number bigger or equal to 2.")
    else
        return 2d^2
    end 
end

"""
    octagonal_color_code_stabilizers(d::Int, color::String)

Return the dictionary for the stabilizers of the periodic octagon color code with distance d.

Notes: We follow the convention in Fig. 2a in https://arxiv.org/pdf/2112.14447.pdf, 
and there should be smarter way of obtaining the dictionary.
"""
function octagonal_color_code_stabilizers(d::Int, color::String)

    if !(color ∈ ("R", "G", "B"))
        error("`color` has to be either `R`, `G` or `B`.")
    end

    if mod(d, 2) != 0 || d < 2
        error("The distance of a periodic octagon color code should be an even number bigger or equal to 2.")
    end    

    dict = Dict{Int64, Vector{Int64}}()
    if color == "R"
        dict[1] = [1, d, 2d^2-d+1, 2d^2] # corner stabilizer
    
        for i = 2 : 2 : (d-2) # (d-2)/2 horiztonal edge stabilizers
            dict[length(dict) + 1] = [i, i+1, i + 2d^2-d, i+1 + 2d^2-d]
        end
        for i = 3d+1 : 4d : (2d^2-5d+1) # (d-2)/2 vertical edge stabilizers
            dict[length(dict) + 1] = [i, i+d-1, i+d, i+2d-1]
        end
    
        # d^2/2 - (d-2) - 1 bulk stabilizers
        for i = d+1 : 4d : (2d^2-3d+1)
            for j = 1 : d/2
                dict[length(dict) + 1] = [i, i+1, i+d, i+d+1] .+ 2(j-1)
            end
    
            if i < (2d^2-3d+1)
                for j = 1 : d/2-1
                    dict[length(dict) + 1] = [i, i+1, i+d, i+d+1] .+ 2(j-1) .+ (2d+1)
                end        
            end
        end
        @assert length(dict) == d^2/2
        
    elseif color == "B"
        # d/2 edge stabilizers
        for i = 1 : 4d : (2d^2-4d+1)
            dict[length(dict)+1] = [i, i+d, i+2d, i+3d, i+d-1, i+2d-1, i+3d-1, i+4d-1]
        end
    
        # d^2/4 - d/2 bulk stabilizers
        for i = 2 : 4d : (2d^2-4d+2)
            for j = 1 : d/2-1
                dict[length(dict)+1] = [i, i+1, i+d, i+d+1, i+2d, i+2d+1, i+3d, i+3d+1] .+ 2(j-1)
            end
        end
        @assert length(dict) == d^2/4
    elseif color == "G"
        # d/2 edge stabilizers 
        for i = 1 : 2 : d-1
            dict[length(dict)+1] = [i, i+1, i+d, i+d+1, i+2d^2-d, i+2d^2-d+1, i+2d^2-2d, i+2d^2-2d+1]
        end
    
        # d^2/4 - d/2 bulk stabilizers
        for i = (1+2d) : 4d : (2d^2-d+1-5d)
            for j = 1 : 2 : d-1
                dict[length(dict)+1] = [i, i+1, i+d, i+d+1, i+2d, i+2d+1, i+3d, i+3d+1] .+ (j-1)
            end
        end
        @assert length(dict) == d^2/4
    end

    return dict
end

"""
    octagonal_color_code_stabilizers(d::Int)

Return the stabilizers for the octagonal color code at distance d
"""
function octagonal_color_code_stabilizers(d::Int)

    dict_R = octagonal_color_code_stabilizers(d, "R")
    dict_G = octagonal_color_code_stabilizers(d, "G")
    dict_B = octagonal_color_code_stabilizers(d, "B")

    # Get rid of the first G and B stabilizers
    pop!(dict_G, 1)
    pop!(dict_B, 1)

    dict = Dict{Int64, Vector{Int64}}()

    for (_, val) in sort(dict_R)
        dict[length(dict)+1] = val
    end

    for (_, val) in sort(dict_G)
        dict[length(dict)+1] = val
    end

    for (_, val) in sort(dict_B)
        dict[length(dict)+1] = val
    end

    return dict        
end


"""
    octagonal_color_code_logicals(d::Int)

Return the logical operator for the octagonal color code at distance d
"""
function octagonal_color_code_logicals(d::Int)
    dict = Dict{Int64, Vector{Int64}}()
    dict[1] = 1 : d
    dict[2] = [1, vcat([[i, i+d] for i in 3d+1 : 4d : (2d^2-5d+1)]...)..., 2d^2-d+1]
    dict[3] = d+1 : 2d
    dict[4] = vcat([[i, i+d] for i in d+1 : 4d : (2d^2-3d+1)]...)
    return dict
end

"""
    octagonal_color_code_normalizers(d::Int)

Return the normalizers (stabilizers and logical operator) for the octagonal color code at distance d
"""
function octagonal_color_code_normalizers(d::Int)
    dict = octagonal_color_code_logicals(d)
    stabs = octagonal_color_code_stabilizers(d)
    for (_, stab) in sort(stabs)
        dict[length(dict)+1] = stab
    end
    return dict
end


"""
    octagonal_color_code_Mq(d::Int)

Return the generator for the periodic octagon color code of distance d in the q subspace
"""
function octagonal_color_code_Mq(d::Int)
    Mq = 2 * diagm(ones(Int64, 2d^2))
    dict = octagonal_color_code_stabilizers(d)
    for (_, value_list) in sort(dict)

        for i in sort(value_list)
            if Mq[i, i] == 2
                for value in value_list
                    Mq[i, value] = 1
                end
                break
            end
        end

#         min_value = minimum(value_list)
    end
    Mq = Mq/√2
    return Mq
end

"""
    octagonal_color_code_Mp(d::Int)

Return the generator for the octagonal color code of distance d in the p subspace
"""
octagonal_color_code_Mp(d) = octagonal_color_code_Mq(d)

"""
    octagonal_color_code_M(d::Int)

Return the generator for the octagonal color code of distance d
"""
function octagonal_color_code_M(d::Int)
    Mq = octagonal_color_code_Mq(d)
    M = BlockDiagonal([Mq, Mq])
    T = basis_transformation(size(Mq)[1])
    return transpose(T) * M * T
end
