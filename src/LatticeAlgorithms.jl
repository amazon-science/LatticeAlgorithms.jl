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

using BlockDiagonals
using LinearAlgebra
using DataStructures
using SparseArrays
using SweepContractor
using PythonCall

const pymatching = PythonCall.pynew()

function __init__()
    # copied from PyPlotCall.jl
    PythonCall.pycopy!(pymatching, pyimport("pymatching"))
    return
end

include("utilities.jl")
export lq_reduce, lll, kz
export islowertriangular, islllreduced, iskzreduced

include("reductions.jl")
export closest_point, shortest_vector, all_closest_points, relevant_vectors
export closest_points

include("lattice_algorithms.jl")
export closest_point_Zn, closest_integer, second_closest_integer
export closest_point_Dn, Dn, Dn_dual, closest_point_Dn_dual
export closest_point_scaled_Zn, closest_point_scaled_Dn, closest_point_scaled_Dn_dual
export An, An_dual, E8, E6, euclidean_dual
export heaviside, next_closest_integer, next_closest_point_Zn!, closest_points_Zn
export closest_points_Dn, closest_points_Dn_dual

include("root_lattices.jl")
export symplectic_dual, bloch_messiah
export general_symplectic, orthogonal_symplectic, get_orthogonal_symplectic_parameters

include("symplectic_utils.jl")
export get_grid_points_2D_lattice, basis_transformation, canonical_form_of_anti_symmetric_matrix
export get_stabilizer_group_from_generators
export get_indicators_from_stabilizers

export Ω_matrix, gram_matrix_of_GKP_lattice_generator, canonize_GKP_lattice_generator
export GKP_logical_operator_generator, GKP_logical_operator_generator_canonical
export distance_X, distance_Y, distance_Z, distance, distances
export gaussian

include("gkp.jl")
export rep_rec, decode_rep_rec
export YY_rep_rec, tlq_YY_rep_rec, decode_YY_rep_rec

include("repetition_codes.jl")

export surface_code_Z_stabilizers, surface_code_X_stabilizers, surface_code_stabilizers
export surface_code_X_logicals, surface_code_Z_logicals
export surface_code_Mq, surface_code_Mp, surface_code_M
export decode_surface_code
export decode_surface_rectangular_GKP_code
export bsv_surface_code
export get_coords_surf_hex, tn_template_surf_hex, tn_surf_hex
export mwms_surface_square, mwms_surface_square_v2
export mwms_surface_hexagonal
include("surface_code.jl")

export rectangular_surface_code_Z_stabilizers
export rectangular_surface_code_X_stabilizers
export rectangular_surface_code_X_logicals
export rectangular_surface_code_Z_logicals
export rectangular_surface_code_Mq
export rectangular_surface_code_Mp
export rectangular_surface_code_M
include("rectangular_surface_code.jl")

export unrotated_surface_code_Z_stabilizers
export unrotated_surface_code_X_stabilizers
export unrotated_surface_code_Z_logicals
export unrotated_surface_code_X_logicals
export unrotated_surface_code_Mq, unrotated_surface_code_Mp
export unrotated_surface_code_M
export decode_unrotated_surface_rectangular_GKP_code
export decode_unrotated_surface_code
export bsv_unrotated_surface_code
export bsv_unrotated_surface_code_qubit
include("unrotated_surface_code.jl")

export decode_concatenated_GKP_code, decode_concatenated_rectangular_GKP_code
include("concatenated_code.jl")

export mwpm
export decoding_graph
export shortest_path, minimum_weight_cycle, mwm, mwms
export shortest_paths
export decoding_graph_v2, mwm_v2, mwms_v2, get_next_matching
include("matching_utils.jl")

export bsv
include("bsv.jl")

export triangular_color_code_num_qubits
export triangular_color_code_stabilizers
export merge_two_stabilizers, triangular_color_code_stabilizers
export triangular_color_code_logicals
export triangular_color_code_normalizers
export triangular_color_code_Mq, triangular_color_code_Mp
export triangular_color_code_M
export get_coords_triangular_color_codes
export tn_template_color_square
export tn_template_color_hex
export tn_color_square, tn_color_hex
export octagonal_color_code_num_qubits
export octagonal_color_code_stabilizers
export octagonal_color_code_stabilizers
export octagonal_color_code_logicals
export octagonal_color_code_normalizers
export octagonal_color_code_Mq, octagonal_color_code_Mp, octagonal_color_code_M
include("color_code.jl")

export coherent_information_pauli_channel, entropy
include("quantum_information_utils.jl")

export get_ws_square, get_ws_non_square
export brute_force_mld_concatenated_square
export brute_force_mld_concatenated_non_square
export sweep_contract_v2!
include("tn_utils.jl")

export gf2_rref, gf2_nullspace, gf2_inverse
export gf2_rref_with_rhs
export gf2_rank, gf2_in_rowspace
export css_stabilizers_from_check_matrix, css_generator_from_check_matrix
export css_logicals, css_distance
export tanner_graph
export cyclic_shift_matrix, circulant_binary_matrix
include("css_ldpc_utils.jl")

export bb_x_matrix, bb_y_matrix
export bb_A_matrix, bb_B_matrix, bb_matrices, bb_check_matrices
export bb_num_data_qubits, bb_num_logicals, bb_parameters
export bb_X_stabilizers, bb_Z_stabilizers, bb_stabilizers
export bb_X_logicals, bb_Z_logicals
export bb_Mq, bb_Mp, bb_M
include("bivariate_bicycle_code.jl")


export hgp_full_rank_repetition_check_matrix, hgp_ring_repetition_check_matrix
export hgp_check_matrices, hgp_num_data_qubits, hgp_num_logicals, hgp_parameters
export hgp_X_stabilizers, hgp_Z_stabilizers, hgp_stabilizers
export hgp_X_logicals, hgp_Z_logicals, hgp_Mq, hgp_Mp, hgp_M
export panteleev_kalachev_c2_parent_check_matrix
include("hyper_graph_product_codes.jl")


export bp_decode
include("bp.jl")

export bp_osd_cs_decode, css_bp_osd_cs_decode
include("bp_osd_cs.jl")

export relay_bp_decode, css_relay_bp_decode
include("bp_relay.jl")

end # module
