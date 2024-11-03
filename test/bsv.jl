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

num_samples = 1e3
for d in [2, 3, 4]
    num_qubits = d^2+(d-1)^2
    stab = unrotated_surface_code_Z_stabilizers(d)
    stabilizers = get_stabilizer_group_from_generators(collect(values(stab)))

    for _ in 1 : num_samples
        ws = 0.6 * randn(num_qubits)
        ws = abs.(ws) # The weight needs to be positive values
        expected_value = sum([prod(ws[stab]) for stab in stabilizers])
        value = bsv(ws)
        if !(10^value ≈ expected_value)
            println(ws)
            println(10^value)
            println(expected_value)
            println()
        end
        @test value ≈ log10(expected_value)
    end
end