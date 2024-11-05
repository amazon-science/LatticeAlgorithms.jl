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

using LinearAlgebra

## Unrotated surface code
@test Set(values(unrotated_surface_code_X_stabilizers(2))) == Set([[1,3,4], [2,3,5]])
@test Set(values(unrotated_surface_code_X_stabilizers(3))) == Set([[2, 4, 5, 7], [8, 10, 13], [7, 9, 10, 12], [6, 9, 11], [3, 5, 8], [1, 4, 6]])
@test Set(values(unrotated_surface_code_Z_stabilizers(2))) == Set([[1,2,3], [3,4,5]])
@test Set(values(unrotated_surface_code_Z_stabilizers(3))) == Set([[4, 6, 7, 9], [10, 12, 13], [5, 7, 8, 10], [2, 3, 5], [9, 11, 12], [1, 2, 4]])
@test Set(values(unrotated_surface_code_X_logicals(2))) == Set([[4, 5]])
@test Set(values(unrotated_surface_code_X_logicals(3))) == Set([[11, 12, 13]])
@test Set(values(unrotated_surface_code_Z_logicals(2))) == Set([[1, 4]])
@test Set(values(unrotated_surface_code_Z_logicals(3))) == Set([[1, 6, 11]])

for d in 4:9
    x_stabilizers = unrotated_surface_code_X_stabilizers(d)
    @test length(x_stabilizers) == length(Set(values(x_stabilizers)))
    @test length(x_stabilizers) == d*(d-1)
    z_stabilizers = unrotated_surface_code_Z_stabilizers(d)
    @test length(z_stabilizers) == length(Set(values(z_stabilizers)))
    @test length(z_stabilizers) == d*(d-1)
end


@test unrotated_surface_code_Mp(2) == [
    1  1  1  0  0
    0  2  0  0  0
    0  0  1  1  1
    0  0  0  2  0
    0  0  0  0  2
]/√2

@test unrotated_surface_code_Mq(2) == [
    1  0  1  1  0
    0  1  1  0  1
    0  0  2  0  0
    0  0  0  2  0
    0  0  0  0  2
]/√2

@test unrotated_surface_code_M(2) == [
    1  0  0  0  1  0  1  0  0  0
    0  1  0  1  0  1  0  0  0  0
    0  0  1  0  1  0  0  0  1  0
    0  0  0  2  0  0  0  0  0  0
    0  0  0  0  2  0  0  0  0  0
    0  0  0  0  0  1  0  1  0  1
    0  0  0  0  0  0  2  0  0  0
    0  0  0  0  0  0  0  2  0  0
    0  0  0  0  0  0  0  0  2  0
    0  0  0  0  0  0  0  0  0  2
]/√2

@test abs(det(unrotated_surface_code_M(2))) ≈ 2
@test abs(det(unrotated_surface_code_M(3))) ≈ 2
@test abs(det(unrotated_surface_code_M(4))) ≈ 2
@test abs(det(unrotated_surface_code_M(5))) ≈ 2



# Test case from BSV https://arxiv.org/pdf/1405.4883.pdf, TABEL I in page 16 and TABLE II in page 17
d = 25
n = d^2 + (d-1)^2

ϵs = 5/100 * ones(n)
prob_I, prob_X = bsv_unrotated_surface_code_qubit(ϵs, zeros(Int, n))


@test (prob_I - log10(1.78283e-27)) < 1e-9
@test (prob_X - log10(5.58438e-57)) < 1e-9

