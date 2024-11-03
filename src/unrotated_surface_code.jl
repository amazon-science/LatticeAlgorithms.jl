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
    unrotated_surface_code_Z_stabilizers(d::Int)

Return the dictionary for the X stabilizers of the unrotated surface code with distance d.

Example:
    >>> unrotated_surface_code_Z_stabilizers(2)
    >>> Dict{Int64, Vector{Int64}} with 2 entries:
        2 => [3, 4, 5]
        1 => [1, 2, 3]

    >>> unrotated_surface_code_Z_stabilizers(3)
    >>> Dict{Int64, Vector{Int64}} with 6 entries:
        5 => [4, 6, 7, 9]
        4 => [10, 12, 13]
        6 => [5, 7, 8, 10]
        2 => [2, 3, 5]
        3 => [9, 11, 12]
        1 => [1, 2, 4]
"""
function unrotated_surface_code_Z_stabilizers(d::Int)
    if d <= 1
        error("`d` should be larger than 1.")
    end    
    Z_dict = Dict{Int64, Vector{Int64}}()
    
    # The stabilizers at the left boundary
    for i in 1 : d-1
        Z_dict[length(Z_dict)+1] = [i, i+1, i+d]
    end
    
    # The stabilizers at the right boundary
    for i in 1 : d-1
        Z_dict[length(Z_dict)+1] = ((2d-1)*(d-1) + i) .+ [1-d, 0 , 1]
    end
    
    # The stabilizers in the bulk
    for i in 1 : d-1
        for j in 1 : d-2
            Z_dict[length(Z_dict)+1] = [d+1, 2d, 2d+1, 3d] .+ (i-1) .+ (2d-1) .* (j-1)
        end
    end
    
    return Z_dict
end

"""
    unrotated_surface_code_X_stabilizers(d::Int)

Return the dictionary for the X stabilizers of the unrotated surface code with distance d.

Example:
    >>> unrotated_surface_code_X_stabilizers(2)
    >>> Dict{Int64, Vector{Int64}} with 2 entries:
        2 => [2, 3, 5]
        1 => [1, 3, 4]

    >>> unrotated_surface_code_X_stabilizers(3)
    >>> Dict{Int64, Vector{Int64}} with 6 entries:
        5 => [2, 4, 5, 7]
        4 => [8, 10, 13]
        6 => [7, 9, 10, 12]
        2 => [6, 9, 11]
        3 => [3, 5, 8]
        1 => [1, 4, 6]
"""
function unrotated_surface_code_X_stabilizers(d::Int)
    if d <= 1
        error("`d` should be larger than 1.")
    end    
    X_dict = Dict{Int64, Vector{Int64}}()
    
    # The stabilizers at the top boundary
    for i in 1 : d-1
        X_dict[length(X_dict)+1] = [1+(2d-1)*(i-1), 1+(2d-1)*(i-1)+d, 1+(2d-1)*(i)]
    end
    
    # The stabilizers at the bottom boundary
    for i in 1 : d-1
        X_dict[length(X_dict)+1] = [d+(2d-1)*(i-1), d+(2d-1)*(i-1)+d-1, d+(2d-1)*(i)]
    end
    
    # The stabilizers in the bulk
    for i in 1 : d-2
        for j in 1 : d-1
            X_dict[length(X_dict)+1] = [
                1+i+(j-1) * (2d-1), 1+i+(j-1) * (2d-1)+d-1,
                1+i+(j-1) * (2d-1)+d, 1+i+(j) * (2d-1)
            ]
        end
    end
    
    return X_dict
end

"""
    unrotated_surface_code_X_logicals(d::Int)

Return the X logical operator for the unrotated surface code with distance d

Example: 
    >>> unrotated_surface_code_X_logicals(3)
    >>> Dict{Int64, Vector{Int64}} with 1 entry:
        1 => [11, 12, 13]
"""
unrotated_surface_code_X_logicals(d::Int) = Dict(1 => collect(d^2+(d-1)^2-d+1:d^2+(d-1)^2))

"""
    unrotated_surface_code_Z_logicals(d::Int)

Return the Z logical operator for the unrotated surface code with distance d

Example: 
    >>> unrotated_surface_code_Z_logicals(3)
    >>> Dict{Int64, Vector{Int64}} with 1 entry:
        1 => [1, 6, 11]
"""
unrotated_surface_code_Z_logicals(d::Int) = Dict(1 => collect(1:(2d-1):d^2+(d-1)^2-d+1))


"""
    unrotated_surface_code_Mq(d::Int)

Return the generator for the surface code of distance d in the q subspace
"""
function unrotated_surface_code_Mq(d::Int)
    Mq = 2 * diagm(ones(Int64, d^2 + (d-1)^2))
    X_dict = unrotated_surface_code_X_stabilizers(d)
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
    unrotated_surface_code_Mp(d::Int)

Return the generator for the surface code of distance d in the p subspace
"""
function unrotated_surface_code_Mp(d::Int)
    Mp = 2 * diagm(ones(Int64, d^2 + (d-1)^2))
    Z_dict = unrotated_surface_code_Z_stabilizers(d)
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
    unrotated_surface_code_M(d::Int)

Return the generator for the surface code of distance d
"""
function unrotated_surface_code_M(d::Int)
    Mq = unrotated_surface_code_Mq(d)
    Mp = unrotated_surface_code_Mp(d)
    M = BlockDiagonal([Mq, Mp])
    T = basis_transformation(size(Mq)[1])
    return transpose(T) * M * T
end


function embed_weights_unrotated(ws::Vector{Float64}, whichtype::String, d::Int)        
    whichtype == "z" && return ws
    GKPind = []
    for i in 1 : d
        push!(GKPind, [i + (2d-1) * (j-1) for j in 1 : d])
        if i < d
            push!(GKPind, [i+d + (2d-1) * (j-1) for j in 1 : d-1])
        end
    end
    GKPind = vcat(GKPind...)
    
    return ws[GKPind]
end    


function get_unrotated_prob(ηs::Vector, whichtype::String, σ::Float64, d::Int; n::Int64=5)
    ws = []
    w2s = []
    for i in 1 : length(ηs)
        cps = closest_points_Zn([(ηs[i] - √π)/2√π], n)
        cp2s = closest_points_Zn([(ηs[i])/2√π], n)

        w = sum([exp(-(ηs[i] - √π - cp[1] * 2√π)^2/(2σ^2)) for cp in cps])
        w2 = sum([exp(-(ηs[i] - cp2[1] * 2√π)^2/(2σ^2)) for cp2 in cp2s])

        push!(ws, w)
        push!(w2s, w2)
    end

    w3s = ws ./ w2s
    
    w3s_embedded = embed_weights_unrotated(w3s, whichtype, d)
    
    Zw3 = bsv(w3s_embedded) + log10(prod(w2s)) - 1/2 * log10(2π * σ^2)

    return Zw3
end

function bsv_unrotated_surface_code(s::Vector, transposeΩMperp::Matrix, σ::Float64; n = 5, subspace="x")
    ηs = transposeΩMperp * s/√(2π) ; 
    d = √(length(s) - 1) / 2 + 1/2
    if !isinteger(d)
        error("The length of ηs has to be 2(d^2+(d-1)^2) for an integer d")
    end
    d = Int(d)    

    function bsv_subspace(η, σ, logical, subspace)

        # For the q subspace, get the prob with no error and with logical error
        prob_I = get_unrotated_prob(η, subspace, σ, d; n=n)
        prob_X = get_unrotated_prob(η - logical, subspace, σ, d; n=n)

        # Determine if there is logical error
        prob_I < prob_X ? lx = 1 : lx = 0

        if lx == 0
            lstar = zeros(d^2 + (d-1)^2)
        else
            lstar = logical
        end
        return lstar
    end

    # Define the candidate shift in the q and p subspaces
    ηs_q = ηs[1:2:end]
    ηs_p = ηs[2:2:end]

    # Define the X and Z logicals
    x_logical = unrotated_surface_code_X_logicals(d)[1]
    z_logical = unrotated_surface_code_Z_logicals(d)[1]

    X = zeros(2(d^2+(d-1)^2))
    X[2x_logical .- 1] .= 1/√2 * √(2π)
    Z = zeros(2(d^2+(d-1)^2))
    Z[2z_logical] .= 1/√2 * √(2π)        

    if subspace == "x"
        lstar = bsv_subspace(ηs_q, σ, X[1:2:end], subspace)
        return -(ηs_q - lstar)
    elseif subspace == "z"
        lstar = bsv_subspace(ηs_p, σ, Z[2:2:end], subspace)
        return -(ηs_p - lstar)
    else
        lstar_q = bsv_subspace(ηs_q, σ, X[1:2:end], "x")
        lstar_p = bsv_subspace(ηs_p, σ, Z[2:2:end], "z")
        rec = ηs
        rec[1:2:end] = rec[1:2:end] - lstar_q
        rec[2:2:end] = rec[2:2:end] - lstar_p
        return -rec
    end
end    

function bsv_unrotated_surface_code_qubit(ϵs::Vector{Float64}, f::Vector{Int64})
        
    #  An aux function for bsv_unrotated_surface_code_qubit(ϵs, f)
    function get_ws_logπf(ϵs::Vector{Float64}, f::Vector{Int64})
        ws = [fe==1 ? (1-ϵse)/ϵse : ϵse/(1-ϵse) for (fe, ϵse) in zip(f, ϵs)]
        logπf = sum([fe==1 ? -log10(wse) : 0 for (fe, wse) in zip(f, ws)])
        return ws, logπf
    end

    N = length(ϵs) # The number of qubits in the unrotated surface code N = d^2 + (d-1)^2

    d = (√(2N-1) + 1)/2 # The distance of the unrotated surface code

    if !isinteger(d)
        error("The length of `ws` and `f` have to be d^2+(d-1)^2 for an integer d")
    end    
    d = round.(Int, d)

    fX = zeros(Int, N)
    fX[1 : (2d-1) : N] .= 1

    ws, logπf = get_ws_logπf(ϵs, f)
    prob_fG = bsv(ws) + logπf + sum(log10.(1 .- ϵs))

    fX2 = mod.(f + fX, 2)
    ws, logπf = get_ws_logπf(ϵs, fX2)
    prob_fXG = bsv(ws) + logπf + sum(log10.(1 .- ϵs))

    return prob_fG, prob_fXG
end