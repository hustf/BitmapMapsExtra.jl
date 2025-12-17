using Test
using BitmapMapsExtras
using BitmapMapsExtras: N0f8, RGBA
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: plot_glyphs!, plot_glyphs, PALETTE_GRGB, RGB
using BitmapMapsExtras: GSTensor, GSTangentBasis, GSVector
using BitmapMapsExtras: Vec2OnGrid, BidirectionOnGrid, 𝐧ₚᵤ!
using BitmapMapsExtras: indices_scattered, indices_on_grid
using BitmapMapsExtras: radial_distance_glyph, pack_glyphs!
using BitmapMapsExtras: indices_on_grid, default_ij_functor, 𝐊!, 𝐊ᵤ!

!@isdefined(hash_image) && include("common.jl")
@testset "Plot tangent basis" begin
    vhash = String[]
    COUNT[] = 0
    gs = GSTangentBasis()
    pts = [CartesianIndex((315, 215))]
    img = plot_glyphs(z_paraboloid(), pts, gs)
    @test is_hash_stored(img, vhash)
    img = background(z_paraboloid())
    pts = indices_on_grid(size(img))
    plot_glyphs!(img, z_paraboloid(), pts, gs)
    @test is_hash_stored(img, vhash)
    img = background(z_sphere_with_bulge())
    plot_glyphs!(img, z_sphere_with_bulge(), pts, gs)

    
end


@testset "Descent vector" begin 
    vhash = ["ae52041380af7ab9d46e34ab9c80a1a8fba62162", "c52e969339d57ceb8a47b89e518612d74e8fa8f8", "f752af27345e6d0f75039e0bab1a444e63c2be07"]
    COUNT[] = 0
    gs = GSVector()
    pts = [CartesianIndex((315, 215))]
    # This creates the default ij functor based on z and gs: With function 𝐧ₚ!.
    img = plot_glyphs(z_paraboloid(), pts, gs)
    @test is_hash_stored(img, vhash)
    img = background(z_paraboloid())
    pts = indices_on_grid(size(img))
    plot_glyphs!(img, z_paraboloid(), pts, gs) 
    @test is_hash_stored(img, vhash)
    img = background(z_sphere_with_bulge(), Δc = 10)
    pts = indices_on_grid(size(img))
    plot_glyphs!(img, z_sphere_with_bulge(), pts, GSVector(; multip = 150, maxg = 150)) 
    @test is_hash_stored(img, vhash)
end

@testset "Descent unit vector" begin 
    vhash = String[]
    COUNT[] = 0
    gs = GSVector()
    pts = [CartesianIndex((315, 215))]
    # 𝐧ₚᵤ! is not the default for a GSVector glyph spec, so we construct 
    # the AbstractIJFunctor here:
    vog = Vec2OnGrid(𝐧ₚᵤ!, z_paraboloid())
    img = plot_glyphs(vog, pts, gs)
    @test is_hash_stored(img, vhash)
    img = background(vog)
    pts = indices_on_grid(vog)
    plot_glyphs!(img, vog, pts, gs)
    @test is_hash_stored(img, vhash)
end


@testset "Curvature" begin
    vhash = String[]
    COUNT[] = 0
    # We don't make a functor object here, because the
    # default for a GSTensor is what we want.
    gs = GSTensor( multip = 12000)
    pts = [CartesianIndex((315, 215))]
    img = plot_glyphs(z_paraboloid(), pts, gs)
    @test is_hash_stored(img, vhash)
    img = background(z_paraboloid())
    pts = indices_on_grid(size(img))
    plot_glyphs!(img, z_paraboloid(), pts, gs)
    @test is_hash_stored(img, vhash)
    img = background(z_sphere_with_bulge())
    pts = indices_on_grid(size(img); Δ = 75)
    plot_glyphs!(img, z_sphere_with_bulge(), pts, GSTensor(; multip = 60000))
    @test is_hash_stored(img, vhash)
    img = background(z_sphere_with_bulge())
    plot_glyphs!(img, z_sphere_with_bulge(), pts, GSTensor(; multip = 60000, direction = 1))
    @test is_hash_stored(img, vhash)
    img = background(z_sphere_with_bulge())
    plot_glyphs!(img, z_sphere_with_bulge(), pts, GSTensor(; multip = 60000, direction = 2))
    @test is_hash_stored(img, vhash)
end

@testset "Unit curvature" begin
    vhash = String[]
    COUNT[] = 0
    gs = GSTensor( multip = 45)
    pts = [CartesianIndex((315, 215))]
    # Normalized length glyphs is not the default BiDirectionOnGrid 
    # given a GSTensor glyph spec.
    # We can supply our own like this: 
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_paraboloid())
    img = plot_glyphs(bdog, pts, gs)
    @test is_hash_stored(img, vhash)
    img = background(bdog)
    pts = indices_on_grid(bdog)
    plot_glyphs!(img, bdog, pts, gs)
    @test is_hash_stored(img, vhash)
end


@testset "Packed" begin
    vhash = String[]
    COUNT[] = 0
    z = z_ellipsoid(;tilt = 0.5)[100:400, 100:400]
    #
    img = background(z)
    pack_glyphs!(img, z,  GSTangentBasis(; halfsize = 15))
    @test is_hash_stored(img, vhash)
    # Packed projected normals
    img = background(z)
    pack_glyphs!(img, z, 
        GSVector(multip = 200, maxg = 200 * 0.2, color = PALETTE_GRGB[4]),
        scatterdist = 4)
    @test is_hash_stored(img, vhash)
    # Packed unit descent
    vog = Vec2OnGrid(𝐧ₚᵤ!, z)
    img = background(vog)
    pack_glyphs!(img, vog, 
        GSVector(multip = 20, maxg = 200 * 0.2, color = PALETTE_GRGB[4]),
        scatterdist = 1)
    @test is_hash_stored(img, vhash)
    # Packed curvature
    img = background(z)
    pack_glyphs!(img, z, 
        GSTensor(multip = 12000, ming = -25, colors=(RGB(0.0,0.049,0.0), RGB(0.0,0.443,1.0))))
    @test is_hash_stored(img, vhash)
    # Packed unit curvature
    bdog = BidirectionOnGrid(𝐊ᵤ!, z)
    img = background(bdog)
    pack_glyphs!(img, bdog, 
        GSTensor(multip = 20, ming = -25, colors=(RGB(0.0,0.049,0.0), RGB(0.0,0.443,1.0))),
        scatterdist = 1)
    @test is_hash_stored(img, vhash)
    # Major curvature direction 
    img = background(z)
    pack_glyphs!(img, z, 
        GSTensor(multip = 12000, ming = -25, direction = 1,
        colors=(RGB(0.0,0.049,0.0), RGB(0.0,0.443,1.0)))) 
    @test is_hash_stored(img, vhash)
    # Minor curvature direction
    img = background(z) 
    pack_glyphs!(img, z, 
        GSTensor(multip = 12000, ming = -25, direction = 2, 
            colors=(RGB(0.0,0.443,1.0), RGB(0.0,0.443,1.0))))
    @test is_hash_stored(img, vhash)
end

@testset "More examples curvature" begin
    vhash = String[]
    COUNT[] = 0
    bdog = BidirectionOnGrid(𝐊!, z_cylinder_offset(π / 6))
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 12000))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊!, z_ellipsoid(; tilt = π / 4))
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 12000, strength = 10))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊!, z_ellipsoid(; tilt = π / 4, a = 0.3))
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 12000, strength = 10))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊!, z_paraboloid())
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 12000))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊!, z_paraboloid(; a = 400, b = 600))
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 10000))
    @test is_hash_stored(img, vhash)
    # 
    bdog = BidirectionOnGrid(𝐊!, z_exp3())
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 7000, maxg = 100, ming=-100))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊!, z_cos(;mult = 50))
    img = 0.3N0f8 .* RGB{N0f8}.background(bdog)).+ RGB{N0f8}(0.08, 0.08, 0.08)
    pack_glyphs!(img, bdog, GSTensor(multip = 17000))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊!, z_ridge_peak_valleys())
    img = 0.3N0f8 .* background(bdog) .+ RGB{N0f8}(0.08, 0.08, 0.08)
    pack_glyphs!(img, bdog, GSTensor(multip = 5000))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊!, z_sphere_with_bulge())
    img = 0.3N0f8 .* background(bdog) .+ RGB{N0f8}(0.08, 0.08, 0.08)
    pack_glyphs!(img, bdog, GSTensor(multip = 35000, direction = 1))
    @test is_hash_stored(img, vhash)
end


@testset "More examples unit curvature" begin
    vhash = String[]
    COUNT[] = 0
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_cylinder_offset(π / 6))
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor())
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_ellipsoid(; tilt = π / 4))
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 30))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_ellipsoid(; tilt = π / 4, a = 0.3))
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 30))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_paraboloid())
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor())
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_paraboloid(; a = 400, b = 600))
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 10000))
    @test is_hash_stored(img, vhash)
    # 
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_exp3())
    img = background(bdog)
    pack_glyphs!(img, bdog, GSTensor(multip = 20, maxg = 100, ming=-100))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_cos(;mult = 50))
    img = 0.3N0f8 .* background(bdog) .+ RGB{N0f8}(0.08, 0.08, 0.08)
    pack_glyphs!(img, bdog, GSTensor(multip = 20))
    @test is_hash_stored(img, vhash)
    #
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_ridge_peak_valleys())
    img = 0.3N0f8 .* background(bdog) .+ RGB{N0f8}(0.08, 0.08, 0.08)
    pack_glyphs!(img, bdog, GSTensor(;multip = 15, strength = 4))
    @test is_hash_stored(img, vhash)
    # 
    bdog = BidirectionOnGrid(𝐊ᵤ!, z_sphere_with_bulge())
    img = 0.3N0f8 .* background(bdog) .+ RGB{N0f8}(0.08, 0.08, 0.08)
    pack_glyphs!(img, bdog, GSTensor(;multip = 15, strength = 4, direction = 2))
    @test is_hash_stored(img, vhash)
    #
    img = 0.3N0f8 .* background(bdog) .+ RGB{N0f8}(0.08, 0.08, 0.08)
    pack_glyphs!(img, bdog, GSTensor(;multip = 15, strength = 4, direction = 1))
    @test is_hash_stored(img, vhash)
end

