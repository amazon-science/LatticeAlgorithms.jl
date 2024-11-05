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
using LinearAlgebra
using Test

"""
    coherent_info_max_mixed_state(px,py,pz)

Expected coherent information for the maximally mixed state.
"""    
function coherent_info_max_mixed_state(px,py,pz)
    p0 = 1 - px - py - pz
    return p0 * log2(p0) + px * log2(px) + py * log2(py) + pz * log2(pz) + 1
end

num_samples = 1e3
for _ in 1 : num_samples
    px, py, pz = rand()/3, rand()/3, rand()/3
    @test coherent_info_max_mixed_state(px,py,pz) ≈ coherent_information_pauli_channel(px,py,pz)
end
