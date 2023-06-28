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

module LatticeAlgorithms

using LatticeAlgorithms
using BlockDiagonals
using LinearAlgebra

include("utilities.jl")

export lq_reduce, lll, kz
export islowertriangular, islllreduced, iskzreduced
include("reductions.jl")

export closest_point, shortest_vector, all_closest_points, relevant_vectors
include("lattice_algorithms.jl")

export closest_point_Zn, closest_points_Zn, closest_integer, second_closest_integer
export closest_point_Dn, Dn, Dn_dual, closest_point_Dn_dual
export closest_point_scaled_Zn, closest_point_scaled_Dn, closest_point_scaled_Dn_dual
export An, An_dual, E8, E6, euclidean_dual
include("root_lattices.jl")


export symplectic_dual, bloch_messiah
export general_symplectic, orthogonal_symplectic, get_orthogonal_symplectic_parameters
include("symplectic_utils.jl")

export get_grid_points_2D_lattice, basis_transformation, canonical_form_of_anti_symmetric_matrix
include("utilities.jl")

export Ω_matrix, gram_matrix_of_GKP_lattice_generator, canonize_GKP_lattice_generator
export GKP_logical_operator_generator, GKP_logical_operator_generator_canonical
export distance_X, distance_Y, distance_Z, distance, distances
include("gkp.jl")

export rep_rec, decode_rep_rec
export YY_rep_rec, tlq_YY_rep_rec, decode_YY_rep_rec
include("repetition_codes.jl")

end # module

