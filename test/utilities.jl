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

using Test
using LatticeAlgorithms

@test basis_transformation(1) == [1 0; 0 1]
@test basis_transformation(2) == [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]

stabilizers = surface_code_Z_stabilizers(3)
stabilizers = get_stabilizer_group_from_generators(collect(values(stabilizers)))

@test Set(stabilizers) == Set([
    [],
    [6, 9],
    [2, 3, 5, 6],
    [9, 2, 3, 5],
    [4, 5, 7, 8],
    [6, 9, 4, 5, 7, 8],
    [2, 3, 6, 4, 7, 8],
    [9, 2, 3, 4, 7, 8],
    [1, 4],
    [6, 9, 1, 4],
    [2, 3, 5, 6, 1, 4],
    [9, 2, 3, 5, 1, 4],
    [5, 7, 8, 1],
    [6, 9, 5, 7, 8, 1],
    [2, 3, 6, 7, 8, 1],
    [9, 2, 3, 7, 8, 1],
])

indicators = get_indicators_from_stabilizers(stabilizers)
@test Set(indicators) == Set([
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, 0, 1, 0, 0, 0, 0],
    [2, 1, 3, 0, 0, 0, 0, 0, 0],
    [2, 1, 1, 0, 1, 0, 0, 0, 0],
    [0, 2, 1, 3, 0, 0, 0, 0, 0],
    [0, 2, 3, 3, 1, 0, 0, 0, 0],
    [2, 3, 2, 3, 0, 0, 0, 0, 0],
    [2, 3, 0, 3, 1, 0, 0, 0, 0],
    [1, 2, 0, 0, 0, 0, 0, 0, 0],
    [1, 2, 2, 0, 1, 0, 0, 0, 0],
    [3, 3, 3, 0, 0, 0, 0, 0, 0],
    [3, 3, 1, 0, 1, 0, 0, 0, 0],
    [1, 0, 1, 3, 0, 0, 0, 0, 0],
    [1, 0, 3, 3, 1, 0, 0, 0, 0],
    [3, 1, 2, 3, 0, 0, 0, 0, 0],
    [3, 1, 0, 3, 1, 0, 0, 0, 0],   
])