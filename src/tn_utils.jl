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
    get_ws_square(η, σ; Nv=5)

Return the weights for the concatenated-square GKP code.

Args:
    η: The candidate errors
    σ: The noise strength
    Nv: The number of closeset points

Returns:
    ws: The vector of weights 
"""
function get_ws_square(η, σ; Nv=5)
    ws = Vector{Float64}()
    for i in 1 : length(η)
        cps = closest_points_Zn([(η[i] - √π)/2√π], Nv)
        w = sum([exp(-(η[i] - √π - cp[1] * 2√π)^2/(2σ^2)) for cp in cps])
        push!(ws, w)
    end
    return ws
end

"""
    brute_force_mld_concatenated_square(ηs_q, σ, logical, full_stabilizers; Nv=5)

Perform brute-force MLD for concatenated square GKP code

Args:
    ηs_q: The candidate errors in the q subspace
    σ: The noise strength
    logical: The logical operator
    full_stabilizers: The full stabilizer group
    Nv: The number of closeset points

Returns:
    lstar: The most probable logical operator
    prob_I: The coset probability for I
    prob_X: The coset probability for X
"""
function brute_force_mld_concatenated_square(ηs_q, σ, logical, full_stabilizers; Nv=5)
    X = zeros(length(ηs_q))
    X[logical] .= 1/√2 * √(2π)
    function get_prob_brute_force_mld_concatenated_square(η)
        ws = get_ws_square(η, σ, Nv=Nv)
        w2s = get_ws_square(η .+ √π, σ, Nv=Nv)
        expected_value = 0
        for stab in full_stabilizers
            prodd = [i ∈ stab ? ws[i] : w2s[i] for i in 1 : length(η)]
            expected_value += prod(prodd)
        end
        return expected_value
    end
    prob_I = get_prob_brute_force_mld_concatenated_square(ηs_q)
    prob_X = get_prob_brute_force_mld_concatenated_square(ηs_q-X)

    prob_I < prob_X ? lstar = X : lstar = zeros(length(X))
    return lstar, prob_I, prob_X
end

"""
    get_ws_non_square(ηs, σ; S_T = [2 0; 1 sqrt(3)] / (12)^(1/4), Nv=5)

Return the weights for the concatenated-square GKP code.

Args:
    ηs: The candidate errors
    σ: The noise strength
    S: The symplectic matrix for the inner one-mode GKP code
    Nv: The number of closeset points

Returns:
    ws: The vector of vector of weights
"""
function get_ws_non_square(ηs, σ; S = [2 1; 0 sqrt(3)] / (12)^(1/4), Nv=5)
    num_qubits = Int(length(ηs)/2)
    ws = []
    S_T = transpose(S)
    X = √π * S * [1, 0]
    Z = √π * S * [0, 1]
    Y = X + Z
    for i in 1 : num_qubits
        x = ηs[2i-1:2i]
        cps_I, cps2_I = closest_points(x    , 2√π * S_T, Nv)
        cps_X, cps2_X = closest_points(x - X, 2√π * S_T, Nv)
        cps_Z, cps2_Z = closest_points(x - Z, 2√π * S_T, Nv)
        cps_Y, cps2_Y = closest_points(x - Y, 2√π * S_T, Nv)

        prob_I = sum([exp(-norm(x-cp)^2/(2σ^2))   for cp in cps_I])
        prob_X = sum([exp(-norm(x-X-cp)^2/(2σ^2)) for cp in cps_X])
        prob_Z = sum([exp(-norm(x-Z-cp)^2/(2σ^2)) for cp in cps_Z])
        prob_Y = sum([exp(-norm(x-Y-cp)^2/(2σ^2)) for cp in cps_Y])
        push!(ws, [prob_I, prob_X, prob_Z, prob_Y])
    end
    return ws
end

"""
    brute_force_mld_concatenated_non_square(ηs, σ, X, Z, full_stabilizers; S = [2 1; 0 √3] / (12)^(1/4), Nv=5)

Perform brute-force MLD for concatenated-non-square GKP code

Args:
    ηs: The candidate errors in the q subspace
    σ: The noise strength
    X: The X logical operator
    Z: The Z logical operator
    full_stabilizers: The full stabilizer group
    S: The symplectic matrix for the inner one-mode GKP code
    Nv: The number of closeset points

Returns:
    lstar: The most probable logical operator
    prob_I: The coset probability for I
    prob_X: The coset probability for X
    prob_Y: The coset probability for Y
    prob_Z: The coset probability for Z
"""
function brute_force_mld_concatenated_non_square(ηs, σ, X, Z, indicators; S = [2 1; 0 √3] / (12)^(1/4), Nv=5)
    function get_prob_brute_force_mld_concatenated_non_square(η)
        ws = get_ws_non_square(η, σ; S = S, Nv=Nv)
        expected_value = 0
        for indicator in indicators
            expected_value += prod([w[i+1] for (w, i) in zip(ws, indicator)])
        end
        return expected_value
    end

    Y = X+Z
    prob_I = get_prob_brute_force_mld_concatenated_non_square(ηs)
    prob_X = get_prob_brute_force_mld_concatenated_non_square(ηs - X)
    prob_Z = get_prob_brute_force_mld_concatenated_non_square(ηs - Z)
    prob_Y = get_prob_brute_force_mld_concatenated_non_square(ηs - Y)

    prob = max(prob_I, prob_X, prob_Y, prob_Z)
    if prob == prob_I
        lstar = zeros(length(ηs))
    elseif prob == prob_X
        lstar = X
    elseif prob == prob_Z
        lstar = Z
    elseif prob == prob_Y
        lstar = Y
    end
 
    return lstar, prob_I, prob_X, prob_Y, prob_Z
end


"""
    sweep_contract_v2!(TN::TensorNetwork, χ::Int, τ::Int; 
start_ind = 1, end_ind = length(TN), MPS_t = [ones(1,1,1)], MPS_i = Int[], report=false)

Modified sweep_contract! from SweepContractor.jl package which allows to output intermediate result

Args:
    TN: The tensor network to be evaluated
    χ: The max bond dim to which truncate the bond dimension
    τ: The max bond dim beyond which performing a truncation
    start_ind: The index of the tensor to start the contraction
    end_ind: The index of the tensor to stop the contraction
    MPS_t, MPS_i: Intermediate results of the tensor network contraction
    report: If reporting the number of truncations
"""

function sweep_contract_v2!(TN::TensorNetwork, χ::Int, τ::Int; 
    start_ind = 1, end_ind = length(TN), MPS_t = [ones(1,1,1)], MPS_i = Int[], report=false)

    # sort!(TN)

    resexp = 0
    count = 0

    for i in start_ind : end_ind
        t = TN[i]
        ind_up = Int[]
        ind_do = Int[]
        for n ∈ t.adj
            if TN[n]>t
                push!(ind_up, n)
            elseif TN[n]<t
                push!(ind_do, n)
            else
                throw(InvalidTNError("Overlapping tensors"))
            end
        end
        sort!(ind_up, by=λ->atan(TN[λ].x-t.x,TN[λ].y-t.y))
        sort!(ind_do, by=λ->atan(TN[λ].x-t.x,t.y-TN[λ].y))
        σ = SweepContractor.permutebetween(t.adj, [ind_do; ind_up])
        t.arr = permutedims(t.arr, σ)
        s = size(t.arr)
        t.arr = reshape(t.arr,(prod(s[1:length(ind_do)]),s[length(ind_do)+1:end]...))

        if isempty(MPS_i)
            MPS_t = SweepContractor.splitMPStensor(MPS_t[1][1]*reshape(t.arr,(size(t.arr)...,1)))
            MPS_i = ind_up
        else
            lo = findfirst(isequal(i), MPS_i)
            hi = findlast(isequal(i), MPS_i)

            isnothing(lo) && throw(InvalidTNError("Disconnected TN"))

            X::Array{Float64} = MPS_t[lo]
            for j ∈ lo+1:hi
                finalsize = (size(X,1),size(X,2)*size(MPS_t[j],2),size(MPS_t[j],3))
                X = reshape(X,(size(X,1)*size(X,2),size(X,3)))*
                    reshape(MPS_t[j],(size(MPS_t[j],1),size(MPS_t[j],2)*size(MPS_t[j],3)))
                X = reshape(X,finalsize)
            end
            X = permutedims(X,[1,3,2])
            M = reshape(t.arr,(size(t.arr,1),prod(size(t.arr)[2:end])))
            X = reshape(
                reshape(X,(size(X,1)*size(X,2),size(X,3)))*M,
                (size(X,1),size(X,2),size(M,2))
            )
            X = permutedims(X,[1,3,2])
            X = reshape(X,(size(X,1),size(t.arr)[2:end]...,size(X,3)))

            MPS_i = [MPS_i[1:lo-1]; ind_up; MPS_i[hi+1:end]]
            if ndims(X)!=2
                MPS_t = [MPS_t[1:lo-1]; SweepContractor.splitMPStensor(X); MPS_t[hi+1:end]]
            elseif isempty(MPS_i)
                MPS_t=[reshape([X[1]],(1,1,1))]
            elseif lo>1
                s = size(MPS_t[lo-1])
                MPS_t[lo-1] = reshape(
                    reshape(MPS_t[lo-1],(s[1]*s[2],s[3]))*X,
                    (s[1],s[2],size(X,2))
                )
                MPS_t = [MPS_t[1:lo-1]; MPS_t[hi+1:end]]
            else
                s = size(MPS_t[hi+1])
                MPS_t[hi+1] = reshape(
                    X*reshape(MPS_t[hi+1],(s[1],s[2]*s[3])),
                    (size(X,1),s[2],s[3])
                )
                MPS_t = [MPS_t[1:lo-1]; MPS_t[hi+1:end]]
            end

            if any(size.(MPS_t,3).>τ)
                count += 1
                SweepContractor.truncMPS!(MPS_t, χ)
                if LinearAlgebra.norm(MPS_t[1])==0
                    return (0.0,typemin(Int))
                end
                h = Int(floor(log2(LinearAlgebra.norm(MPS_t[1]))))
                resexp += h
                MPS_t[1] /= exp2(h)
            end
        end
    end

    report && println("Number of truncations: $count")

    if end_ind == length(TN)
        res = MPS_t[1][1];
        if res == 0.0
            return (0.0, typemin(Int64));
        end
        h = Int(floor(log2(abs(res))));
        return res/exp2(h),resexp+h, count
    else
        return MPS_t, MPS_i, resexp, count
    end
end