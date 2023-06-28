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

@testset "lattice_algorithms" begin
    include("lattice_algorithms.jl")
end


@testset "reductions" begin
    include("reductions.jl")
end


@testset "root_lattices" begin
    include("root_lattices.jl")
end


@testset "symplectic_utils" begin
    include("symplectic_utils.jl")
end


@testset "utilities" begin
    include("utilities.jl")
end

@testset "gkp" begin
    include("gkp.jl")
end


@testset "repetition_codes" begin
    include("repetition_codes.jl")
end
