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
    closest_integer(x::Float64)

Return the closest integers to x. In case of a tie, choose the integer with the smallest absolute value.
"""
closest_integer(x::Float64) = mod(x, 1)==0.5 ? trunc(Int, x) : round(Int, x) 

"""
    function second_closest_integer(x::Float64)

Return the second closest integers to x. In case of a tie, choose the integer with the smallest absolute value.

Ref: http://neilsloane.com/doc/Me83.pdf
"""
function second_closest_integer(x::Float64)
    if x==0.0
        return 1
    elseif mod(abs(x),1) <= 0.5
        return round(Int,sign(x)) * (abs(closest_integer(x))+1)
    else
        return round(Int,sign(x)) * (abs(closest_integer(x))-1)
    end
end

# Zn

"""
    closest_point_Zn(x::Vector)

Return the closest point for the given vector x in the Zn lattice.

Note that here the generator of the Zn lattice is taken to be the identity matrix

"""
closest_point_Zn(x::Vector) = closest_integer.(x)

"""
    closest_point_scaled_Zn(x::Vector, scale::Number)

Return the closest point for the given vector x in the scaled Zn lattice.

Note that here the generator of the scaled Zn lattice is taken to be a diagonal matrix
"""
closest_point_scaled_Zn(x::Vector, scale::Number) = closest_point_Zn(x./scale) .* scale

"""
    closest_point_scaled_Zn(x::Vector, scales::Vector)

Return the closest point for the given vector x in the Zn lattice where each axis has different scales.

Note that here the generator of the scaled Zn lattice is taken to be a diagonal matrix
"""
closest_point_scaled_Zn(x::Vector, scales::Vector) = closest_point_Zn(x./scales) .* scales


# Dn

"""
    Dn(N::Int)

Return a representation for the Dn lattice.
"""
function Dn(N::Int)
    @assert N>=1
    if N==1
        return Matrix(2I, 1, 1)
    end
    
    M = zeros(N, N)
    for i = 1 : N-1
        M[i, i], M[i, i+1] = 1, -1
    end
    M[N, N-1], M[N, N] = 1, 1
    
    return M
end


"""
    closest_point_scaled_Dn(x::Vector; scales::Vector)

Return the closest points for the given x in the Dn lattice where each axis has different scale

Note that here the generator of the D_n lattice is the one returned by Dn(N::Int) function scaled by `scales`
"""
function closest_point_scaled_Dn(x::Vector, scales::Vector)
    fx = closest_integer.(x./scales)

    if mod(sum(fx), 2) == 0
        return fx.*scales
    else
        wx = second_closest_integer.(x./scales)
        δx1 = (x .- fx .* scales).^2
        δx2 = (x .- wx .* scales).^2
        diff = abs.(δx1 - δx2)
        ind = findfirst(y -> y == minimum(diff), diff)
        gx = copy(fx)
        gx[ind] = wx[ind]
        return gx.*scales
    end
end


"""
    closest_point_Dn(x::Vector)

Return the closest point for the given x in the Dn lattice

Note that here the generator of the D_n lattice is the one returned by Dn(N::Int) function
"""
closest_point_Dn(x::Vector) = closest_point_scaled_Dn(x, ones(length(x)))


"""
    closest_point_scaled_Dn(x::Vector, scale::Number)

Return the closest point for the given x in the scaled Dn lattice
    
Note that here the generator of the D_n lattice is the one returned by Dn(N::Int) function scaled by `scale`
"""
closest_point_scaled_Dn(x::Vector, scale::Number) = closest_point_scaled_Dn(x, scale * ones(length(x)))


"""
    euclidean_dual(M::Matrix)

Return an Euclidean dual of the given generator matrix 
"""
euclidean_dual(M::Matrix) = Matrix(inv(transpose(M)))


# Dn_dual 

"""
    Dn_dual(N::Int)

Return a representation for the Euclidean dual of NxN Dn lattice

"""
Dn_dual(N::Int) = euclidean_dual(Dn(N))

"""
    closest_point_Dn_dual(x::Vector, scales::Vector)

Return the closest point for the given x in the dual of Dn lattice

Note that here the generator of the dual of the D_n lattice is the one returned by Dn_dual(N::Int) function scaled by `scales`
"""
function closest_point_scaled_Dn_dual(x::Vector, scales::Vector)
    n = length(x)
    y = x ./ scales
    r1 = 0.5 * ones(n)
    y0 = closest_integer.(y) .* scales
    y1 = (closest_integer.(y - r1) + r1) .* scales
    d0 = norm(x - y0)
    d1 = norm(x - y1)
    if d0 < d1
        return y0
    else
        return y1
    end
end


"""
    closest_point_Dn_dual(x::Vector)

Return the closest point for the given x in the dual of Dn lattice

Note that here the generator of the dual of the D_n lattice is the one returned by Dn_dual(N::Int) function
"""
closest_point_Dn_dual(x::Vector) = closest_point_scaled_Dn_dual(x, ones(length(x)))

"""
    closest_point_Dn_dual(x::Vector, scale::Number)

Return the closest point for the given x in the dual of Dn lattice

Note that here the generator of the dual of the D_n lattice is the one returned by Dn_dual(N::Int) function scaled by `scale`
"""
closest_point_scaled_Dn_dual(x::Vector, scale::Number) = closest_point_scaled_Dn_dual(x, scale * ones(length(x)))

# A-type

"""
    An(N::Int)

Return a representation for the Dn lattice.

Note: By definition, the Euclidean Gram matrix of the A-type root lattice, namely An(N)*transpose(An(N)), is equal to Matrix(1I, N, N) + ones(N, N). 
See Eq. 53 in Chapter 4 of Conway-Sloane
"""
An(N::Int) = Matrix(1I, N, N) - (√(N+1)+1)/N * ones(N, N)

"""
    An_dual(N::Int)

Return a representation for the Euclidean dual of NxN An lattice
"""
An_dual(N::Int) = euclidean_dual(An(N))


# E-type

"""
    E8()

Return a representation for the E8 lattice.

Ref: Eq. 99 in Chapter 4 of Conway-Sloane. 
"""
E8() = [
    2 0 0 0 0 0 0 0;
    -1 1 0 0 0 0 0 0;
    0 -1 1 0 0 0 0 0;
    0 0 -1 1 0 0 0 0;
    0 0 0 -1 1 0 0 0;
    0 0 0 0 -1 1 0 0;
    0 0 0 0 0 -1 1 0;
    1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/2;
]


"""
    E6()

Return a representation for the E6 lattice.

Ref: Eq. 7a in http://neilsloane.com/doc/Me108.pdf
"""
E6() = [
    0 √3 0 0 0 0;
    0 0 0 √3 0 0;
    1 0 1 0  1 0;
    -3/2 -√3/2 0 0 0 0;
    0 0 -3/2 -√3/2 0 0;
    -1/2 √3/2 -1/2 √3/2 -1/2 √3/2
]
