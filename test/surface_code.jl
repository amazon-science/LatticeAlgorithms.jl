using LatticeAlgorithms
using LinearAlgebra

@test surface_code_Mp(3) == [
    1  0  0  1  0  0  0  0  0
    0  1  1  0  1  1  0  0  0
    0  0  2  0  0  0  0  0  0
    0  0  0  1  1  0  1  1  0
    0  0  0  0  2  0  0  0  0
    0  0  0  0  0  1  0  0  1
    0  0  0  0  0  0  2  0  0
    0  0  0  0  0  0  0  2  0
    0  0  0  0  0  0  0  0  2
]/√2

@test surface_code_Mq(3) == [
    1  1  0  1  1  0  0  0  0
    0  1  1  0  0  0  0  0  0
    0  0  2  0  0  0  0  0  0
    0  0  0  2  0  0  0  0  0
    0  0  0  0  1  1  0  1  1
    0  0  0  0  0  2  0  0  0
    0  0  0  0  0  0  1  1  0
    0  0  0  0  0  0  0  2  0
    0  0  0  0  0  0  0  0  2
]/√2

@test surface_code_M(3) == [
    1  0  1  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0
    0  1  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0
    0  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0
    0  0  0  1  0  1  0  0  0  1  0  1  0  0  0  0  0  0
    0  0  0  0  2  0  0  0  0  0  0  0  0  0  0  0  0  0
    0  0  0  0  0  2  0  0  0  0  0  0  0  0  0  0  0  0
    0  0  0  0  0  0  2  0  0  0  0  0  0  0  0  0  0  0
    0  0  0  0  0  0  0  1  0  1  0  0  0  1  0  1  0  0
    0  0  0  0  0  0  0  0  1  0  1  0  0  0  1  0  1  0
    0  0  0  0  0  0  0  0  0  2  0  0  0  0  0  0  0  0
    0  0  0  0  0  0  0  0  0  0  2  0  0  0  0  0  0  0
    0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  1
    0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0
    0  0  0  0  0  0  0  0  0  0  0  0  0  2  0  0  0  0
    0  0  0  0  0  0  0  0  0  0  0  0  0  0  2  0  0  0
    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2  0  0
    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2  0
    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2    
]/√2

@test abs(det(surface_code_M(3))) ≈ 2
@test abs(det(surface_code_M(5))) ≈ 2

@test surface_code_X_logicals(3)[1] == [1,4,7]
@test surface_code_Z_logicals(3)[1] == [1,2,3]

@test Set(values(surface_code_Z_stabilizers(3))) == Set([[2,3,5,6], [4,5,7,8], [1,4], [6,9]])
@test Set(values(surface_code_X_stabilizers(3))) == Set([[1,2,4,5], [5,6,8,9], [2,3], [7,8]])
@test Set(values(surface_code_stabilizers(3))) == Set([[2, 8],[9, 11, 15, 17], [4, 6, 10, 12], [8, 10, 14, 16], [13, 15], [12, 18], [3, 5], [1, 3, 7, 9]])

@test Set(values(surface_code_Z_stabilizers(5))) == Set(
    [
        [2,3,7,8], [4,5,9,10], [6,7,11,12], [8,9,13,14], [1,6], [11, 16],
        [12,13,17,18], [14,15,19,20], [16,17,21,22], [18,19,23,24], [10, 15], [20, 25]
    ]
)

@test Set(values(surface_code_X_stabilizers(5))) == Set(
    [
        [1,2,6,7], [3,4,8,9], [7,8,12,13], [9,10,14,15], [2,3], [4,5],
        [11,12,16,17], [13,14,18,19], [17,18,22,23], [19,20,24,25], [21,22], [23,24]
    ]
)

function test_surface_code_stabilizers(d)
    
    # Test surface_code_stabilizers via generator matrix
    M1 = surface_code_M(d)

    M2 = 2 * diagm(ones(Int64, 2d^2))
    dict = surface_code_stabilizers(d)
    for (_, value_list) in dict
        min_value = minimum(value_list)
        for value in value_list
            M2[min_value, value] = 1
        end
    end
    M2 = M2/√2

    @test M1 ≈ M2
end

for d in 3:2:19
    test_surface_code_stabilizers(d)
end
