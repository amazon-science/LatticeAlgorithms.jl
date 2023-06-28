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
using LinearAlgebra
using BlockDiagonals
using LatticeAlgorithms

@test_throws ErrorException symplectic_dual(Matrix(I, 3, 3))
@test basis_transformation(1) == [1 0; 0 1]
@test basis_transformation(2) == [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]

## massive tests

function test_orthogonal_symplectic(N::Integer, samples::Integer=Int(1e2))
    ω = [0 1; -1 0]
    Ω = BlockDiagonal([ω for _ in 1:N])

    S = orthogonal_symplectic(N)

    elapsed_time = @elapsed for _ in 1 : samples
        args = rand(1,N^2)
        S0 = S(args)

        @test det(S0) ≈ 1
        @test S0 * Ω * transpose(S0) ≈ Ω
        @test S0 * transpose(S0) ≈ Matrix(1.0I,2*N,2*N)
        @test transpose(S0) * S0 ≈ Matrix(1.0I,2*N,2*N)
    end
    # println("test_orthogonal_symplectic N=$N, elapsed_time=$elapsed_time")
end

for N in 1 : 6
    test_orthogonal_symplectic(N)
end


function test_general_symplectic(N::Integer, method::String, samples::Integer=Int(1e2))
    ω = [0 1; -1 0]
    Ω = BlockDiagonal([ω for _ in 1:N])

    S = general_symplectic(N, method = method)

    elapsed_time = @elapsed for _ in 1 : samples
        args = rand(1,2*N^2 + N)
        S0 = S(args)

        @test det(S0) ≈ 1
        @test S0 * Ω * transpose(S0) ≈ Ω
    end
    # println("test_general_symplectic_$method N=$N, elapsed_time=$elapsed_time")
end
for N in 1 : 6
    test_general_symplectic(N, "generator")
    test_general_symplectic(N, "Bloch_Messiah")
end


function test_bloch_messiah(N::Integer, samples::Integer = Int(1e2))
    S = general_symplectic(N); 
    Ω = BlockDiagonal([[0 1; -1 0] for _ in 1 : N])

    elapsed_time = @elapsed for _ in 1 : samples
        args = rand(2N^2 + N) ; 
        S0 = S(args); 
    
        O1, Z, O2 = bloch_messiah(S0)
        
        # check that O1 * Z * O2 = S0
        @test O1 * Z * O2 ≈ S0
    
        # check that transpose(O1) * O1 = 1    
        @test transpose(O1) * O1 ≈ Matrix(1I, 2N, 2N)
    
        # check that O1 * transpose(O1) = 1
        @test O1 * transpose(O1) ≈ Matrix(1I, 2N, 2N)
        
        # check that transpose(O2) * O2 = 1    
        @test transpose(O2) * O2 ≈ Matrix(1I, 2N, 2N)
    
        # check that O2 * transpose(O2) = 1
        @test O2 * transpose(O2) ≈ Matrix(1I, 2N, 2N)

        # check that O1 is symplectic
        @test O1 * Ω * transpose(O1) ≈ Ω
    
        # check that O2 is symplectic
        @test O2 * Ω * transpose(O2) ≈ Ω        

    end
    # println("test_bloch_messiah N=$N", elapsed_time=$elapsed_time")

end

for N in 1 : 9
    test_bloch_messiah(N)
end