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

using LatticeAlgorithms
using Test
using LinearAlgebra

M = Dn(4)

@test Ω_matrix(M) == [0 1 0 0; -1 0 0 0; 0 0 0 1; 0 0 -1 0]

@test gram_matrix_of_GKP_lattice_generator(M) == [0 1 0 0; -1 0 1 -1; 0 -1 0 2; 0 1 -2 0]

@test canonize_GKP_lattice_generator(M)[2] == [0 2 0 0; -2 0 0 0; 0 0 0 1; 0 0 -1 0]

@test GKP_logical_operator_generator(M) == [1 -1 0 0; 0 1 0 0; 1/2 -1/2 1/2 -1/2; -1/2 1/2 1/2 1/2]

@test GKP_logical_operator_generator_canonical(M) == [0 1 0 0; -1/2 -1/2 1/2 1/2; 0 0 1 1; 0 1 -1 0]

@test [distance_X(M), distance_Y(M), distance_Z(M)] ≈ [2.506628274631, 2.506628274631, 2.506628274631]
@test distances(M) ≈ [2.506628274631, 2.506628274631, 2.506628274631]
@test distance(M) ≈ 2.506628274631

@test canonical_form_of_anti_symmetric_matrix(gram_matrix_of_GKP_lattice_generator(M))[2] == [0 2 0 0; -2 0 0 0; 0 0 0 1; 0 0 -1 0]
@test abs(det(canonical_form_of_anti_symmetric_matrix(gram_matrix_of_GKP_lattice_generator(M))[1])) == 1

@test_throws ErrorException Ω_matrix(Dn(3))

@test_throws ErrorException gram_matrix_of_GKP_lattice_generator([1 1; 2 -2.1])

@test_throws ErrorException canonical_form_of_anti_symmetric_matrix([1 1; 2 -2.1])