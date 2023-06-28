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

struct TLQBasis 
    T:: Matrix{Int64}
    L::Matrix
    Q::Matrix
end


"""
    function lq_reduce(M::Matrix)

Return L, Q such that M = L * Q. 

Here L is a lower triangular matrix with positive diagonal elements, and Q is an orthogonal matrix.
"""
function lq_reduce(M:: Matrix)
    Q, R = qr(transpose(M)) 
    Q, L = Matrix(transpose(Q)), Matrix(transpose(R))

    # The above produce M = L * Q where L is lower triangular. 
    # To make sure L has positive diagonals, we do the following
    for i in 1 : size(M)[1]
        if L[i, i] < 0
            L[:, i] = -L[:, i] # flip the column of L
            Q[i, :] = -Q[i, :] # flip the row of Q
        end
    end
    return TLQBasis(Matrix(1I, size(M)), L, Q)
end


"""
    lll(M::Matrix; δ::Float64=3/4)

Implement the LLL reduction to return T, L, Q such that M = T * L * Q 

Here T is a unimodular matrix, Q is an orthogonal matrix and L is a lower triangular matrix which is the LLL reduced basis.

Ref1: https://github.com/christianpeel/LLLplus.jl/blob/master/src/lll.jl
Ref2: https://www.ant.uni-bremen.de/sixcms/media.php/102/10740/SPM_2011_Wuebben.pdf, pg 10. This is copied from their MATLAB code attached. 

Note that we assume M is a square matrix 
"""
function lll(M::Matrix; δ::Float64=3/4)
    if !(0.25 < δ < 1.0)
        error("δ must be between 1/4 and 1.");
    end

    N = size(M)[1]
    B = copy(transpose(M)) # The column vectors of B are basis vectors
    Q, R = qr(B) # B=QR, where Q is an orthogonal matrix, R is upper triangular 

    Q = Matrix(Q) 

    T = Matrix(1I, N, N) # unimodular NxN matrix

    l = 2

    while l <= N
        for k = l-1 : -1 : 1
            mu = round(R[k,l]/R[k,k]) # abs(R[k,l])>0.5*abs(R[k,k])
            if abs(mu) > 0
                B[:, l] = B[:, l] - mu * B[:, k]
                R[1:k, l] = R[1:k, l] - mu * R[1:k, k]
                T[:, l] = T[:, l] - mu * T[:, k]
            end
        end


        # Lovasz condition 
        len = norm(R[l-1:l,l]) 
        if δ * abs(R[l-1,l-1])^2 > len^2

            # swapping of columns l-1 and l in B, T and R
            B[:,[l-1 l]]   = B[:,[l l-1]];
            T[:,[l-1 l]]   = T[:,[l l-1]];
            R[1:l,[l-1 l]] = R[1:l,[l l-1]];

            # reconstruction of upper triangular structure by Givens rotation 
            # mutliplication with matrix Theta achieves R(l,l-1) = 0
            c = R[l-1,l-1]/len # len = ||R[l-1:l,l-1]|| after swapping
            s = R[l,l-1]/len
            Theta = [[c' s]; [-s c]]

            R[l-1:l,l-1:end] = Theta * R[l-1:l,l-1:end];
            Q[:,l-1:l] = Q[:,l-1:l] * Theta' ;
            l = max(l-1,2);
        else
            l = l+1;
        end

    end

    B = Matrix(transpose(B)) # The row vectors of B are basis vectors 
    T = Matrix(transpose(T)) # Such that B = T*M
    Q = Matrix(Q)
    R = Matrix(R)
    
    for i in 1 : size(R)[1]
        if R[i, i] < 0
            Q[:, i] = -Q[:, i] # flip the column of Q
            R[i, :] = -R[i, :] # flip the row of R
        end
    end

    T = inv(T)
    L = Matrix(transpose(R)) 
    Q = Matrix(transpose(Q))

    @assert T * L * Q ≈ M
    @assert islowertriangular(L)
    @assert Q * transpose(Q) ≈ Matrix(1I, size(Q))
    @assert transpose(Q) * Q ≈ Matrix(1I, size(Q))
    @assert round.(Int, T) ≈ T
    @assert abs(det(T)) ≈ 1

    return TLQBasis(round.(Int, T), L, Q)

end

"""
    kz(M::Matrix)

Implement the KZ reduction to return T, L, Q such that M = T * L * Q

Here T is a unimodular matrix, Q is an orthogonal matrix and L is a lower triangular matrix which is the KZ reduced basis.

Ref: http://www.cas.mcmaster.ca/~qiao/publications/ZQW11.pdf

The KZ reduction is essentially solving N shortest vector problems, recursively. 

"""
function kz(M::Matrix; δ::Float64=3/4)
    B = copy(transpose(M))
    N = size(B)[1]

    Q, R = qr(B)
    Q = Matrix(Q) 

    Z = Matrix(1I, N,N)
    for k=1:N-1
        lll_basis = lll(Matrix(transpose(R[k:end, k:end])), δ=δ)
        T1, L1, Q1 = lll_basis.T, lll_basis.L, lll_basis.Q

        z = round.(Int, transpose(inv(T1)) * shortest_vector(L1))

        for j=N-k+1:-1:2
            if z[j] != 0
                d, a, b = gcd_extended(z[j-1], z[j])
                T2 = [z[j-1]/d -b;
                      z[j]/d    a]
                
                z[j-1] = d
                R[1:j+k-1, j+k-2:j+k-1]= R[1:j+k-1, j+k-2:j+k-1]*T2
                
                # eliminate R[j+k-1, j+k-2]
                len2 = sqrt(sum(R[j+k-2:j+k-1, j+k-2].^2))
                c = R[j+k-2, j+k-2]/len2
                s = R[j+k-1, j+k-2]/len2
                Theta = [[c' s]; [-s c]]
                
                R[j+k-2:j+k-1,j+k-2:end] = Theta*R[j+k-2:j+k-1,j+k-2:end]
                Z[:,j+k-2:j+k-1] = Z[:,j+k-2:j+k-1]*T2
                # Q[:,j+k-2:j+k-1] = Q[:,j+k-2:j+k-1] * Theta' ;
            end
        end
    end
    # size reduce
    for l=2:N
        # reduce l-th column of B
        for k=l-1:-1:1
            mu = round(Int, R[k,l]/R[k,k])
            if abs(mu)>0
                R[1:k,l] .-= mu .* R[1:k,k]
                Z[:,l] .-= mu .* Z[:,k]
            end
        end
    end

    Z = round.(Int,Z)
    T = Matrix(transpose(inv(Z)))
    L = Matrix(LowerTriangular(transpose(R)))
    Q = Matrix(inv(L) * transpose(Z) * M)

    # To make sure L has positive diagonals
    for i in 1 : size(L)[1]
        if L[i, i] < 0
            Q[i, :] = -Q[i, :]
            L[:, i] = -L[:, i]
        end
    end

    # @assert T * L * Q ≈ M
    # @assert islowertriangular(L)
    # @assert Q * transpose(Q) ≈ Matrix(1I, size(Q))
    # @assert transpose(Q) * Q ≈ Matrix(1I, size(Q))
    # @assert round.(Int, T) ≈ T
    # @assert abs(det(T)) ≈ 1
    
    return TLQBasis(round.(Int, T), L, Q)
end


"""
    islowertriangular(M::Matrix)

Return True if the input matrix M is lower triangular, otherwise return False.
"""
islowertriangular(M::Matrix) = (Matrix(LowerTriangular(M)) ≈ M)


"""
    function islllreduced(M::Matrix)

Return True if the input matrix M is LLL reduced, otherwise return False.

For the condition of being LLL reduced, see Eq. 4-9 in the following Ref
https://publications.lib.chalmers.se/records/fulltext/14990/local_14990.pdf

"""
function islllreduced(M::Matrix)
    @assert islowertriangular(M) == true

    N = size(M)[1]

    for i = 1 : N
        if i == N
            return true
        end

        if norm(M[i, i:end]) > 2/√3 * norm(M[i+1, i:end])
            println("The 1st LLL condition failed")
            return false
        end
        
        for k = i+1 : N
            if abs(M[k,i]) > abs(M[i,i])/2 
                println("The 2nd LLL condition failed at k=$k")
                return false
            end
        end
    end
end


"""
    iskzreduced(L::Matrix, verbose::Bool=false)

Return True if the input matrix L is KZ reduced, otherwise return False.
If verbose is true, then the failed KZ condition will be printed out.

For the condition of being KZ reduced, see Eq. 4-7 in the following Ref
https://publications.lib.chalmers.se/records/fulltext/14990/local_14990.pdf

"""
function iskzreduced(L::Matrix, verbose::Bool=false)
    @assert islowertriangular(L) == true

    N = size(L)[1]

    for i = 1 : N
        if i == N
            return true
        end

        sv = shortest_vector(L[i:end, i:end])
        if !(norm(L[i, i]) ≈ norm(sv))
            if verbose
                println("The 1st KZ condition failed at L[$i,$i]")
            end
            return false
        end
        
        for k = i+1 : N
            if !(abs(L[k,i]) <= abs(L[i,i])/2)
                if verbose
                    println("The 2nd KZ condition failed as abs(L[$k,$i]) > abs(L[$i,$i])/2")
                end
                return false
            end
        end
    end

end

