using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: paint_convexity_rank!, indices_on_grid, RGB
using BitmapMapsExtras: convexity_rank, BidirectionOnGrid, 𝐊!
using BitmapMapsExtras: is_bidirec_vect_positive, N0f8

!@isdefined(is_hash_stored) && include("common.jl")

@testset "Convexity rank" begin
    r = TestMatrices.r
    flatval = 0.0003
    # 1    Convex - convex           Green
    bog = BidirectionOnGrid(𝐊!, z_paraboloid(; a= 0.6r, b = 0.5r))
    @test convexity_rank(bog(400, 400), flatval) == 1
    # 2     Convex - flat            Bleak green
    bog = BidirectionOnGrid(𝐊!, -z_cylinder(0.1))
    @test convexity_rank(bog(400, 400), flatval) == 2
    # 3     Convex - concave         Blue
    bog = BidirectionOnGrid(𝐊!, z_paraboloid(;a = 0.6r, b = -0.4r))
    @test convexity_rank(bog(400, 400), flatval) == 3
    # 4     Flat - flat              Grey
    bog = BidirectionOnGrid(𝐊!, z_plane())
    @test convexity_rank(bog(400, 400), flatval) == 4
    # 5    Flat - concave            Bleak red
    bog = BidirectionOnGrid(𝐊!, z_cylinder(0.1))
    @test convexity_rank(bog(400, 400), flatval) == 5
    # 6    Concave - concave         Red
    bog = BidirectionOnGrid(𝐊!, z_ellipsoid())
    @test convexity_rank(bog(400, 400), 0.0002) == 6
end


@testset "Edge case, curvature sign" begin
    i, j = 500, 232
    bog = BidirectionOnGrid(𝐊!, z_ellipsoid())
    K = bog(i,j)
    @test !is_bidirec_vect_positive(K[:, 1])
    @test !is_bidirec_vect_positive(K[:, 2])
end

@testset "Curvature types" begin
    vhash = ["143b2a749b1e721b1257fe690e7f7649c487412b", "a744293b958ae236e793bdd3a3a4ff22a4479029", "b06d1171f773ac0a73af5d8d65740b78ffe7512d", "ae99f490315aab9e8ac51aef2265b70007620655", "2d5d12b79f1705286dcae14f2bd5c00b7e221f10", "8242248a1d96deba6f886d2a5eb9c70f0ca96b56"]
    COUNT[] = 0
    # We compromise on this value
    flatval = 0.0003
    pts = TestMatrices.R
    r = TestMatrices.r
    # 1    Convex - convex           Green
    img =  paint_convexity(z_paraboloid(; a= 0.6r, b = 0.5r); flatval)
    @test is_hash_stored(img, vhash)
    # 2     Convex - flat            Bleak green
    img =  paint_convexity(-z_cylinder(0.1); flatval)
    @test is_hash_stored(img, vhash)
    # 3     Convex - concave         Blue
    img =  paint_convexity(z_paraboloid(;a= 0.6r, b = -0.4r); flatval)
    @test is_hash_stored(img, vhash)
    # 4     Flat - flat              Grey
    img =  paint_convexity(z_plane(); flatval)
    @test is_hash_stored(img, vhash)
    # 5    Flat - concave            Bleak red
    img =  paint_convexity(z_cylinder(0.1); flatval)
    @test is_hash_stored(img, vhash)
    # 6    Concave - concave         Red
    img =  paint_convexity(z_ellipsoid(); flatval)
    @test is_hash_stored(img, vhash)
end

