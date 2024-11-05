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

############ rectangular (rotated) surface code
"""
    rectangular_surface_code_Z_stabilizers(dx::Int, dz::Int)

Return the dictionary for the Z stabilizers of the rectangular surface code with distances dx and dz.

Example: 
    >>> rectangular_surface_code_Z_stabilizers(3, 5)
    >>> Dict{Int64, Vector{Int64}} with 6 entries:
          5 => [8, 9, 13, 14]
          4 => [6, 7, 11, 12]
          6 => [10, 15]
          2 => [2, 3, 7, 8]
          3 => [4, 5, 9, 10]
          1 => [1, 6]
"""
function rectangular_surface_code_Z_stabilizers(dx::Int, dz::Int)
    if mod(dx, 2) != 1 || mod(dz, 2) != 1
        error("The distances of the surface code should be odd numbers greater than 1.")
    end
    
    dx == dz && return surface_code_Z_stabilizers(dx)
    
    Z_dict = Dict{Int64, Vector{Int64}}()
    num_cols = Int((dz + 1) / 2)
    num_rows = dx - 1
    for row_idx in 1:2:num_rows
        for col_idx in 1:1
            stab_idx = col_idx + (row_idx - 1) * num_cols
            qubit_idx = 2 * col_idx - 2 + (row_idx - 1) * dz
            Z_dict = merge!(Z_dict, 
                Dict(stab_idx => [qubit_idx + 1, qubit_idx + dz + 1])
            )
        end
        for col_idx in 2:num_cols
            stab_idx = col_idx + (row_idx - 1) * num_cols
            qubit_idx = 2 * col_idx - 2 + (row_idx - 1) * dz 
            Z_dict = merge!(Z_dict, 
                Dict(stab_idx => [qubit_idx, qubit_idx + 1, qubit_idx + dz, qubit_idx + dz + 1])
            )
        end
    end
    for row_idx in 2:2:num_rows
        for col_idx in 1:num_cols-1
            stab_idx = col_idx + (row_idx - 1) * num_cols
            qubit_idx = 2 * col_idx - 1 + (row_idx - 1) * dz 
            Z_dict = merge!(Z_dict, 
                Dict(stab_idx => [qubit_idx, qubit_idx + 1, qubit_idx + dz, qubit_idx + dz + 1])
            )
        end
        for col_idx in num_cols:num_cols
            stab_idx = col_idx + (row_idx - 1) * num_cols
            qubit_idx = 2 * col_idx - 1 + (row_idx - 1) * dz
            Z_dict = merge!(Z_dict, 
                Dict(stab_idx => [qubit_idx, qubit_idx + dz])
            )
        end
    end
    return Z_dict
end

"""
    rectangular_surface_code_X_stabilizers(dx::Int, dz::Int)

Return the dictionary for the X stabilizers of the rectangular surface code with distances dx and dz.

Example: 
    >>> rectangular_surface_code_X_stabilizers(3, 5)
    >>> Dict{Int64, Vector{Int64}} with 8 entries:
          5 => [3, 4, 8, 9]
          4 => [7, 8, 12, 13]
          6 => [13, 14]
          7 => [4, 5]
          2 => [11, 12]
          8 => [9, 10, 14, 15]
          3 => [2, 3]
          1 => [1, 2, 6, 7]
"""
function rectangular_surface_code_X_stabilizers(dx::Int, dz::Int)
    if mod(dx, 2) != 1 || mod(dz, 2) != 1
        error("The distances of the surface code should be odd numbers greater than 1.")
    end
    
    dx == dz && return surface_code_X_stabilizers(dx)
    
    X_dict = Dict{Int64, Vector{Int64}}()
    num_cols = dz - 1
    num_rows = Int((dx + 1) / 2)
    for col_idx in 1:2:num_cols
        for row_idx in 1:num_rows-1
            stab_idx = row_idx + (col_idx - 1) * num_rows
            qubit_idx = (row_idx - 1) * 2 * dz + col_idx
            X_dict = merge!(X_dict, 
                Dict(stab_idx => [qubit_idx, qubit_idx + 1, qubit_idx + dz, qubit_idx + dz + 1])
            )
        end
        for row_idx in num_rows:num_rows
            stab_idx = row_idx + (col_idx - 1) * num_rows
            qubit_idx = (row_idx - 1) * 2 * dz + col_idx
            X_dict = merge!(X_dict, 
                Dict(stab_idx => [qubit_idx, qubit_idx + 1])
            )
        end
    end
    for col_idx in 2:2:num_cols
        for row_idx in 1:1
            stab_idx = row_idx + (col_idx - 1) * num_rows
            qubit_idx = (row_idx - 1) * 2 * dz - dz + col_idx
            X_dict = merge!(X_dict, 
                Dict(stab_idx => [qubit_idx + dz, qubit_idx + dz + 1])
            )
        end
        for row_idx in 2:num_rows
            stab_idx = row_idx + (col_idx - 1) * num_rows
            qubit_idx = (row_idx - 1) * 2 * dz - dz + col_idx
            X_dict = merge!(X_dict, 
                Dict(stab_idx => [qubit_idx, qubit_idx + 1, qubit_idx + dz, qubit_idx + dz + 1])
            )
        end
    end
    return X_dict
end

"""
    rectangular_surface_code_Z_logicals(dx::Int, dz::Int)

Return the Z logical operator for the rectangular surface code with distances dx and dz

Example: 
    >>> rectangular_surface_code_Z_logicals(3, 5)
    >>> Dict{Int64, Vector{Int64}} with 1 entry:
        1 => [1, 2, 3, 4, 5]
"""
rectangular_surface_code_Z_logicals(dx::Int, dz::Int) = Dict(1 => collect(1:dz))

"""
    rectangular_surface_code_X_logicals(dx::Int, dz::Int)

Return the X logical operator for the rectangular surface code with distances dx and dz

Example: 
    >>> rectangular_surface_code_X_logicals(3, 5)
    >>> Dict{Int64, Vector{Int64}} with 1 entry:
        1 => [1, 6, 11]
"""
rectangular_surface_code_X_logicals(dx::Int, dz::Int) = Dict(1 => collect(1:dz:dx*dz))
    
    
"""
    rectangular_surface_code_Mq(dx::Int, dz::Int)

Return the generator for the rectangular surface code of distance dx and dz in the q subspace
"""
function rectangular_surface_code_Mq(dx::Int, dz::Int)
    Mq = 2 * diagm(ones(Int64, dx * dz))
    X_dict = rectangular_surface_code_X_stabilizers(dx, dz)
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
    rectangular_surface_code_Mp(dx::Int, dz::Int)

Return the generator for the rectangular surface code of distance dx and dz in the p subspace
"""
function rectangular_surface_code_Mp(dx::Int, dz::Int)
    Mp = 2 * diagm(ones(Int64, dx * dz))
    Z_dict = rectangular_surface_code_Z_stabilizers(dx, dz)
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
    rectangular_surface_code_M(dx::Int, dz::Int)

Return the generator for the rectangular surface code of distance dx and dz
"""
function rectangular_surface_code_M(dx::Int, dz::Int)
    Mq = rectangular_surface_code_Mq(dx, dz)
    Mp = rectangular_surface_code_Mp(dx, dz)
    M = BlockDiagonal([Mq, Mp])
    T = basis_transformation(size(Mq)[1])
    return transpose(T) * M * T
end
