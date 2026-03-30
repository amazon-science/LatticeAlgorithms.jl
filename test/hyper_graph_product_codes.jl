using Test
using LinearAlgebra

using LatticeAlgorithms

@testset "hyper_graph_product_codes" begin
    @testset "toric code as HGP of ring repetition" begin
        # Standard HGP benchmark: toric code from repetition parents.
        d = 3
        H = hgp_ring_repetition_check_matrix(d)
        HX, HZ = hgp_check_matrices(H, H)

        @test size(H) == (d, d)
        @test size(HX, 2) == 2 * d^2
        @test size(HZ, 2) == 2 * d^2

        # CSS commutation condition.
        @test all(iszero, mod.(HX * transpose(HZ), 2))

        # For d = 3, the toric code is [[18, 2, 3]].
        params = hgp_parameters(H, H; compute_distance=true)
        @test params.n == 18
        @test params.k == 2
        @test params.d == 3
    end

    @testset "generic cyclic parent helper" begin
        H = circulant_binary_matrix(7, [0, 1, 3])

        @test size(H) == (7, 7)
        @test all(sum(H, dims=1) .== 3)
        @test all(sum(H, dims=2) .== 3)

        # Check circulant structure: each row is a cyclic shift of the first.
        for r in 0:6
            @test H[r + 1, :] == circshift(H[1, :], r)
        end
    end

    @testset "Panteleev-Kalachev C2 parent" begin
        H = panteleev_kalachev_c2_parent_check_matrix()

        @test size(H) == (31, 31)
        @test all(sum(H, dims=1) .== 3)
        @test all(sum(H, dims=2) .== 3)

        # The first row matches h(x) = 1 + x^2 + x^5.
        @test H[1, :] == [1, 0, 1, 0, 0, 1, zeros(Int64, 25)...]

        # Check circulant structure.
        for r in 0:30
            @test H[r + 1, :] == circshift(H[1, :], r)
        end
    end

    @testset "Panteleev-Kalachev C2 HGP parameters" begin
        H = panteleev_kalachev_c2_parent_check_matrix()
        HX, HZ = hgp_check_matrices(H, H)

        # CSS commutation condition.
        @test all(iszero, mod.(HX * transpose(HZ), 2))

        # Literature parameters: [[1922, 50, 16]].
        # In unit tests we check n and k only; exact d=16 is too expensive for CI.
        params = hgp_parameters(H, H; compute_distance=false)
        @test params.n == 1922
        @test params.k == 50

        # For this self-product construction, both check matrices are (3,6)-regular in the paper.
        @test all(sum(HX, dims=2) .== 6)
        @test all(sum(HZ, dims=2) .== 6)
    end
end