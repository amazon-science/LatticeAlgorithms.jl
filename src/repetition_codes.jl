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
    rep_rec(N::Int)

Return the generator for the generalized Tesseract code with N modes.

The N-mode generalized Tesseract code is defined to be the concatenation 
of N-mode repetition code (with XX stabilizers here) and the rectangular 
GKP code with ratio N^(1/2).

"""
function rep_rec(N::Int)
    M = 2 * Matrix(1I, 2N, 2N)    
    for i = 1 : N-1
        M[2i-1, 2i-1:2i+2] = [1 0 1 0]
    end
    M = M/√2 * BlockDiagonal([[N^(1/4) 0; 0 N^(-1/4)] for _ in 1:N])
    return M
end

"""
    decode_rep_rec(ξ::Vector)

Return the closest point for the given syndrome in the symplectic dual lattice of the generalized Tesseract code.

Here we find the closest point for the given ξ in the lattice generated by √(2π) * Mperp, where Mperp is the
generator for the symplectic dual lattice in the canonical basis. The corresponding logical operator can be 
obtained as 
    u = inv(transpose(√(2π) * Mperp)) * y
    u = mod.(round.(Int, u[1:2]), 2)
"""
function decode_rep_rec(ξ::Vector)
    dim = length(ξ)
    @assert mod(dim, 2)==0
    
    N = div(dim,2) # number of modes

    λ_q = √(2π) * N^(1/4)/√2 
    λ_p = √(2π) * √2/N^(1/4) 
    
    ξ_q = ξ[1:2:2N-1]
    ξ_p = ξ[2:2:2N]
    
    y_q = closest_point_Zn(ξ_q/λ_q) * λ_q
        
    y_p = closest_point_Dn_dual(ξ_p/λ_p) * λ_p

    y = zeros(2N)
    y[1:2:end], y[2:2:end] = y_q, y_p

    return y
end


"""
    YY_rep_rec(N::Int)

Return the generator for the generalized Tesseract code with 2N modes.

The 2N-mode generalized Tesseract code is defined to be the concatenation 
of YY stabilizer with two copies of N-mode generalized Tesseract codes.

"""
function YY_rep_rec(N::Int)
    M = Matrix(
        BlockDiagonal(
            [N^(1/4)/√2 * Dn(N), 
             N^(1/4)/√2 * Dn(N),
             √2/N^(1/4) * Matrix(1I, N, N), 
             √2/N^(1/4) * Matrix(1I, N, N)
            ]
        )
    )
    
    # r is the YY stabilizer
    r = vcat([0 for _ in 1:2N], [N^(-1/4) for _ in 1:2N])
    r[1] = r[N+1] = N^(1/4)
    r = 1/√2 * r
    M[end, :] = r
    
    T2 = basis_transformation(2N)
    M = transpose(T2) * M * T2
    
    return M  
end


"""
    tlq_YY_rep_rec(N::Int)

Return T, L, Q, r where T * L * Q = Mperp is the logical operator for the 2N-mode generalized Tesseract code, 
    and r is the glueing operator for the code. 

"""
function tlq_YY_rep_rec(N::Int)

    function Dn2(N)
        M = Matrix(0I, N, N)
        for i = 1 : N-1
            M[i, i:i+1] = [1 1]
        end
        M[end, end] = 2
        return M
    end
    
    function Dn_dual2(N)
        M = Matrix(1.0I, N, N)
        M[end, :] .= 1/2
        return M
    end

    M = YY_rep_rec(N)
    Mperp = GKP_logical_operator_generator_canonical(M); 
    
    # The easy-to-decode logical operator
    Mperp2 = Matrix(BlockDiagonal([N^(1/4)/√2 * Dn2(2N), N^(-1/4)*√2 * Dn_dual2(2N)]))
    
    r = [[0 for _ in 1:2N-1]..., N^(1/4)/√2, [N^(-1/4)/√2 for _ in 1:N]..., [0 for _ in 1:N]...]
    
    Mperp2[3N, :] = r
    
    L = Mperp2

    T2 = basis_transformation(2N)    
    Mperp3 = transpose(T2) * Mperp2 * T2    
    R = Mperp * inv(Mperp3)
    

    T = R * transpose(T2)
    Q = T2
    
    @assert T * L * Q ≈ Mperp
    @assert abs(det(T)) ≈ 1
    @assert round.(Int, T) ≈ T
    @assert Q * transpose(Q) ≈ Matrix(1I, 4N, 4N)
    @assert transpose(Q) * Q ≈ Matrix(1I, 4N, 4N)
    
    return T, L, Q, r
end

"""
    decode_YY_rep_rec(ξ, T2, L2, Q2, r2) 

Return the closest point for the given point with respect to the 2N-mode generalized Tesseract code in linear time.

Note `_, _, Q2, r2 = tlq_YY_rep_rec(N)` where N = Int(length(ξ)/4).
"""
function decode_YY_rep_rec(ξ, Q2, r2) 
    ξ2 = Q2 * ξ/√(2π)
    N = Int(length(ξ)/4)

    rs = [[0 for _ in 1:4N], r2]
    ds = []
    zs = []
    λ_q, λ_p = N^(1/4)/√2, √2/N^(1/4)
    for r in rs
        x = ξ2-r
        x_q, x_p = x[1:2N], x[2N+1:4N]
        z_q = closest_point_Dn(x_q/λ_q)
        z_p = closest_point_Dn_dual(x_p/λ_p)
        z = vcat(z_q * λ_q, z_p * λ_p)

        push!(ds, norm(z - x))
        push!(zs, z)
    end    
    
    ind_min = findmin(ds)[2] # index for the vector with smallest distance

    z = zs[ind_min] + rs[ind_min] # shift the origin of the vector with the smallest distance

    y = √(2π) * transpose(z) * Q2
    return vec(y)
end