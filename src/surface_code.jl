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