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
    closest_point(x::Vector, M::Union{Matrix, TLQBasis})

Return the closest point for the given vector x in the lattice generated by M.

For the cloesest point of x, denoted as y, one could extract its coordinates in the lattice via 
    round.(Int, transpose(y) * inv(M))
where the rounding is important since the coordinates are always integer.

Ref 1: https://publications.lib.chalmers.se/records/fulltext/14990/local_14990.pdf
"""
closest_point(x::Vector, M::Matrix) = closest_point(x, lll(M))

function closest_point(x::Vector, M::TLQBasis)
    T, L, Q = M.T, M.L, M.Q
    @assert size(T) == (length(x), length(x))

    u_hat_temp = decode(L, x=vec(transpose(x) * transpose(Q)), problem="closest_point")
    u_hat = transpose(u_hat_temp) * inv(T)
    # return u_hat
    return vec(u_hat * T*L*Q)
end



"""
    shortest_vector(M::Matrix)

Return the shortest vector in the lattice generated by M.

Ref 1: https://publications.lib.chalmers.se/records/fulltext/14990/local_14990.pdf
"""
shortest_vector(M::Matrix) = shortest_vector(lll(M))

function shortest_vector(M::TLQBasis)
    T, L, Q = M.T, M.L, M.Q

    u_hat_temp = decode(L, problem="shortest_vector")
    u_hat = transpose(u_hat_temp) * inv(T)
    return vec(u_hat * (T*L*Q))
end


"""
    decode(L:: Matrix; x::Union{Vector, Nothing}=nothing, problem::String = "closest_point")

Decode the lattice generator L for various lattice problems.

Valid `problem` are:
    - "closest_point": return u_hat such that u_hat * L is the cloest lattice point to the given point x. 
    - "shortest_vector": return u_hat such that u_hat * L is the shortest vector of the lattice
    - "all_closest_points": return U_hat ≡ {u_hat∈Z^N: u_hat * L is the cloest lattice point to the given point x}
    - "kissing_number": 

Here L is an NxN lower triangular matrix with positive element, which is the generator matrix of a lattice. x is the point of interest which is an N dimensional real-valued vector.
"""
function decode(L:: Matrix; x::Union{Vector, Nothing}=nothing, problem::String = "closest_point")
    L = inv(L)
    sgn(z) = z<=0 ? -1 : 1 # if z<=0, sgn(z)=-1, else sgn(z)=1

    if !islowertriangular(L)
        error("L has to be lower triangular.")
    end

    function update!(k, dir, u, offset, e)
        """For a given layer with dimension k, and the direction dir , update the examined lattice point u, 
        the offset to the given layer offset, and the distance to the given layer y, 
        """
        if dir=="up" # move up to layers with one more dimension
            u[k] = u[k] + offset[k]
            y = (e[k][k] - u[k]) / L[k, k]
            offset[k] = -offset[k] - sgn(offset[k])
        else # dir=="down", move down to layers with one less dimension
            u[k] = round(Int, e[k][k])
            y = (e[k][k] - u[k]) / L[k, k]
            offset[k] = sgn(y)
        end
        return u, y, offset
    end


    N = size(L)[1] # dimension of A
    bestdist = Inf # current distance record
    k = N # dimension of the examined layer (a subspace of the lattice)

    u = [NaN for _ in 1:N] # candicate lattice point, for not explored coordinates, set to NaN
    offset = [NaN for _ in 1:N] # offset to next layers, for those not explored, set to NaN
    dist = [NaN for _ in 1:N] # distance to the layers, for those not explored, set to NaN
    dist[k] = 0

    e = Dict() # used to compute hat_un

    if problem == "closest_point"
        if x == nothing
            error("For closest_point problem, the input vector x is needed.")
        end
        e[k] = vec(transpose(x) * L) 
        ϵ = 0.0 # threshold to determine the shortest distance
    elseif problem == "shortest_vector"
        e[k] = zeros(1,N)
        ϵ = 0.0
    elseif problem == "all_closest_points"
        if x == nothing
            error("For all_closest_points problem, the input vector x is needed.")
        end
        e[k] = vec(transpose(x) * L) 
        ϵ = 1e-10
        global U_hat = [] # set of all closest_points
    elseif problem == "kissing_number"
        
    end

    u, y, offset = update!(k, "down", u, offset, e)

    while true
        newdist = dist[k] + y^2

        if newdist < (1+ϵ) * bestdist 
            if k != 1 # case A
                if k-1 ∉ keys(e)
                    e[k-1] = [] # just initialize the key k-1
                end
                e[k-1] = [e[k][i] - y * L[k, i] for i in 1:k-1] 
                k = k-1
                dist[k] = newdist 
                u, y, offset = update!(k, "down", u, offset, e)
            else # case B
                if problem == "closest_point" || (problem == "shortest_vector" && newdist > 0.0) # exclude 0 as a solution for shortest_vector problem
                    global u_hat = copy(u)  # best lattice point so far
                    bestdist = newdist
                    k = k+1 
                elseif problem == "all_closest_points"
                    push!(U_hat, copy(u))
                    bestdist = min(newdist, bestdist)
                end
                u, y, offset = update!(k, "up", u, offset, e)
            end
        else # case C
            if k==N
                if problem == "closest_point" || problem == "shortest_vector"
                    return u_hat
                elseif problem == "all_closest_points"
                    return U_hat
                end
            else
                k = k+1
                u, y, offset = update!(k, "up", u, offset, e)
            end
        end
    end
end

"""
    all_closest_points(x::Vector, M::Union{Matrix, TLQBasis})

Return the all the closest points for the given vector x in the lattice generated by M.

Ref 1: https://publications.lib.chalmers.se/records/fulltext/14990/local_14990.pdf
"""
all_closest_points(x::Vector, M::Matrix) = all_closest_points(x, lll(M))

function all_closest_points(x::Vector, M::TLQBasis)

    T, L, Q = M.T, M.L, M.Q
    @assert size(T) == (length(x), length(x))

    U_hat_temp = decode(L, x=vec(transpose(x) * transpose(Q)), problem="all_closest_points")
    U_hat = [transpose(u_hat_temp) * inv(T) for u_hat_temp in U_hat_temp]
    return [vec(u_hat * (T*L*Q)) for u_hat in U_hat]
end



"""
    function relevant_vectors(M::Matrix)

Return Ntil, which is the set of Voronoi-relevant vectors N(Λ) for the lattice generated by M. 
Note that The Voronoi cell at 0∈Λ is defined as Vor_M(n=0)≡{x∈R^N: |x|<|x-c| for ∀c∈N(Λ)}
"""
function relevant_vectors(M::Matrix)    

    # Get Mtil for M
    N = size(M)[1]
    Mtil = []
    for i = 1:2^N-1
       z = bitstring(i)
       z = [1//2 * parse(Int, zz) for zz in collect(z[end-N+1:end])]
    
       push!(Mtil, vec(transpose(z)*M))
    end

    Ntil = []
    for s in Mtil
        H_hat = all_closest_points(s, M)
        if length(H_hat)==2
            push!(Ntil, 2 * (H_hat[1]-s))
        end
    end

    return Ntil 
end