using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: paint_steepness_rank!, steepness_radian
using BitmapMapsExtras: N0f8, indices_on_grid, RGB, mark_at!

!@isdefined(is_hash_stored) && include("common.jl")

@testset "Steepness radian" begin
    @test isapprox(steepness_radian([i for i = 1:5, j = 1:5]), π / 4, atol = 1e-5)
    @test isapprox(steepness_radian([j for i = 1:5, j = 1:5]), π / 4, atol = 1e-5)
    @test isapprox(steepness_radian([-i for i = 1:5, j = 1:5]), π / 4, atol = 1e-5)
    @test isapprox(steepness_radian([9 for i = 1:5, j = 1:5]), 0, atol = 1e-15)
    @test isapprox(steepness_radian([j for i = 1:5, j = 1:5]), π / 4, atol = 1e-5)
    @test isapprox(steepness_radian([i / √3 for i = 1:5, j = 1:5]), π / 6, atol = 1e-5)
    # At cos(π / 3), slope is -π / 6. Pixels make inexact.
    pt = CartesianIndex(499, 499 + 250)
    @test isapprox(steepness_radian(z_sphere(), pt), π / 6, atol = 1e-2)
end

@testset "Paint 22.5° steepness interval" begin
    vhash = String[]
    COUNT[] = 0
    img = background(z_sphere())
    paint_steepness_rank!(img, z_sphere())
    @test is_hash_stored(img, vhash)
    img = background(z_wavy())
    paint_steepness_rank!(img, z_wavy())
    @test is_hash_stored(img, vhash)
end


@testset "Paint 30° steepness interval" begin
    vhash = ["9ff916ac03f7c01565a3353812f818a29144454d", "5101b17c3df5b1f6816ad2a43099a28fa3c254ba", "470128a6176577a94e7f79178913587df140dd23", "57bc257803e50640b4564e80c9add9af0ec1ffc5", "4c3bc79dad254c2ffd157b0393560d0ce239cf5b", "b6e54c144b2796d91475dd52c1b3f02b76450909", "427f89be372ed33b7e63011eaac2ea7671654a2a", "b05fe3dbdbf7e4f3a4dc9aaf83eb4299da485cd0", "a7fb3482e059211e7bfd5540d4cfe10b3bfc86c1", "705398144ad554f8a6ca03ad1f7928c217db0ca8", "98c35961a4338a531e877eadaea3066dd2ff8b43"]
    COUNT[] = 0
    vzf = [z_cos,
        () -> z_cylinder(π / 6),
        () -> z_cylinder_offset(π / 3),
        z_ellipsoid,
        z_exp3,
        z_paraboloid,
        z_ridge_peak_valleys,
        z_sphere,
        z_plane,
        z_wavy, 
        z_sphere_with_bulge]
    for fz in vzf
        z = fz()
        img = background(z)
        paint_steepness_rank!(img, z; steepness_interval_rad = π / 6)
        @test is_hash_stored(img, vhash)
    end
end