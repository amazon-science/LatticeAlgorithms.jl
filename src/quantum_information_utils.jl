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
function coherent_information_pauli_channel(px, py, pz, ρ) 

Return the coherence information for a Pauli channel with respect 
    to an initial state ρ. px, py, pz are the parameters of the 
    Pauli channel. 
"""
function coherent_information_pauli_channel(
    px::Number, 
    py::Number, 
    pz::Number;
    ρ::Matrix = Matrix(I,2,2)/2
) 

    @assert size(ρ) == (2,2)

    X = [0 1; 1 0]
    Y = [0 -1im; 1im 0]
    Z = [1 0; 0 -1]

    p0 = 1 - px - py - pz
    Np = p0 * ρ + px * X * ρ * X + py * Y * ρ * Y + pz * Z * ρ * Z

    Npc = [[p0 √(px*p0)*tr(X*ρ) √(py*p0)*tr(Y*ρ) √(pz*p0)*tr(Z*ρ)]; 
           [√(px*p0)*tr(X*ρ) px -1im*√(px*py)*tr(Z*ρ) 1im*√(pz*px)*tr(Y*ρ)]; 
           [√(py*p0)*tr(Y*ρ) 1im*√(px*py)*tr(Z*ρ) py -1im*√(py*pz)*tr(X*ρ)]; 
           [√(pz*p0)*tr(Z*ρ) -1im*√(pz*px)*tr(Y*ρ) 1im*√(py*pz)*tr(X*ρ) pz]]

    Ic = entropy(Np) - entropy(Npc)
    return Ic
end

"""
function entropy(ρ::Matrix)

Return the von Neumann entropy for a given density matrix ρ.
"""
function entropy(ρ::Matrix)
    eigvalρ = eigvals(ρ)
    eigvalρ = eigvalρ[eigvalρ.>0]
    return -sum(eigvalρ .* log2.(eigvalρ))
end
