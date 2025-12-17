# Test glyph packing
using Test
using BitmapMapsExtras
using BitmapMapsExtras: PALETTE_GRGB, plot_glyph_given_value!, apply_color_by_coverage!
using BitmapMapsExtras: radial_distance_glyph, GSVector, norm, plot_glyphs!, plot_glyphs
using BitmapMapsExtras: placements_and_values, pack_glyphs!, Vec2OnGrid, BidirectionOnGrid
using BitmapMapsExtras: plot_glyphs_given_values!, 𝐧ₚ!, 𝐊!, coarse_radius_for_plotting
using BitmapMapsExtras: GSTensor
using BitmapMapsExtras: coverage_fitting_image, MMatrix
using BitmapMapsExtras: mark_at!, N0f8
using BitmapMapsExtras.TestMatrices
!@isdefined(is_hash_stored) && include("common.jl")

@testset begin "Mark circumference of a vector glyph"
    vhash = String[]
    COUNT[] = 0
    img = fill(RGB{N0f8}(0.624,0.667,0.663), 200, 200)
    cov = [0f0 for i = 1:200, j=1:200]
    gs = GSVector(multip = 1, color = PALETTE_GRGB[3], maxg = 100)
    ptexist = CartesianIndex(100, 100)
    v = [80.0, -20.0] #(x, y)
    plot_glyph_given_value!(cov, ptexist, v, gs)
    apply_color_by_coverage!(img, cov, gs)
    mark_at!(img, ptexist; side = 1)
    rc = coarse_radius_for_plotting(gs, v)
    for α = 0:0.02:2π
        # Point on rough outer circle
        Δi = -Int(round(rc * sin(α)))
        Δj = Int(round(rc * cos(α)))
        npt = ptexist + CartesianIndex((Δi, Δj))
        mark_at!(img, npt; side = 1)
        # 
        # From pt to ptexist
        u = ptexist.I .- npt.I
        # The angle for the ray is calculated as if we're in a 'y up' c.s.
        αu = atan(-u[1], u[2]) 
        # αu ≈ α + π
        rptexist = radial_distance_glyph(gs, v, αu + π)
        ptint = ptexist - CartesianIndex(Int.(round.(u .* rptexist ./ norm(u))))
        mark_at!(img, ptint)
    end
    @test is_hash_stored(img, vhash)
end

@testset begin "Vector glyph crashless placements and values"
    ppts = [CartesianIndex(100, 100), 
            CartesianIndex(100, 110)]
    gs = GSVector(multip = 1, maxg = 100)
    function b!(v, x)
        v .= [80.0, -20.0]
        v
    end
    d = Vec2OnGrid(b!, zeros(Float64, 199, 199))
    # Visual check
    img = fill(RGB{N0f8}(0.624,0.667,0.663), 200, 200)
    plot_glyphs!(img, d, ppts, gs)
    display_if_vscode(img)
    # Only the first glyph is accepted (basic overlap)
    @test placements_and_values(d, gs, ppts) == ([ppts[1]], [[80.0, -20.0]])
    #
    ppts = [CartesianIndex(100, 100), 
            CartesianIndex(100, 130)]
    img = fill(RGB{N0f8}(0.624,0.667,0.663), 200, 200)
    plot_glyphs_given_values!(img, ppts, [[80.0, -20.0], [80.0, -20.0]], gs)
    # Only the first glyph is accepted (glancing overlap)
    @test length(placements_and_values(d, gs, ppts)[1]) == 1
    #
    ppts = [CartesianIndex(100, 100), 
            CartesianIndex(100, 140)]
    # Visual check
    img = fill(RGB{N0f8}(0.624,0.667,0.663), 200, 200)
    plot_glyphs!(img, d, ppts, gs)
    display_if_vscode(img)
    # Both glyphs are accepted (but the order may well be another than in ppts)
    @test length(placements_and_values(d, gs, ppts)[1]) == 2
    @test length(placements_and_values(d, gs, ppts)[1]) == 2
    @test placements_and_values(d, gs, ppts)[2][1] == [80, -20]
end

@testset begin "pack vector glyphs no collision"
    vhash = String[]
    COUNT[] = 0
    img = fill(PALETTE_GRGB[1], 200, 200)
    b(x, y) = [80.0, -20.0]
    @test_throws MethodError pack_glyphs!(img, b, GSVector())
    d = Vec2OnGrid(𝐧ₚ!, z_paraboloid()[400:599, 400:599])
    @test d(3, 99) == [-0.0037305447382605097, -0.36559338434952804]
    gs = GSVector(multip = 30 / norm(d(9,99)), maxg = 100, color = PALETTE_GRGB[2])
    pack_glyphs!(img, d, gs)
    @test is_hash_stored(img, vhash)
end


@testset begin "Mark circumference of a bidirectional vector glyph"
    vhash = String[]
    COUNT[] = 0
    img = fill(RGB{N0f8}(0.624,0.667,0.663), 200, 200)
    gs = GSTensor(multip = 0.8, direction = 1, maxg = 200, ming = -200)
    cov = coverage_fitting_image(img, gs)
    @test cov isa Matrix{Float32}
    ptexist = CartesianIndex(100, 100)
    v = [80.0, 20.0] #(x, y)
    K = MMatrix{2, 2, Float64}(hcat(v, v))
    plot_glyph_given_value!(cov, ptexist, K, gs)
    apply_color_by_coverage!(img, cov, gs)
    mark_at!(img, ptexist; side = 1)
    rc = coarse_radius_for_plotting(gs, K)
    for α = 0:0.005:2π
        # Point on rough outer circle
        Δi = -Int(round(rc * sin(α)))
        Δj = Int(round(rc * cos(α)))
        npt = ptexist + CartesianIndex((Δi, Δj))
        mark_at!(img, npt; side = 1)
        # 
        # From pt to ptexist
        u = ptexist.I .- npt.I
        # The angle for the ray is calculated as if we're in a 'y up' c.s.
        αu = atan(-u[1], u[2]) 
        # αu ≈ α + π
        rptexist = radial_distance_glyph(gs, v, αu + π)
        ptint = ptexist - CartesianIndex(Int.(round.(u .* rptexist ./ norm(u))))
        mark_at!(img, ptint)
    end
    img
    @test is_hash_stored(img, vhash)
end

@testset begin "Mark circumference of a tensor glyph"
    vhash = String[]
    COUNT[] = 0
    img = fill(RGB{N0f8}(0.624,0.667,0.663), 500, 500)
    gs = GSTensor(multip = 2, direction = 1:2, maxg = 500, ming = -500)
    cov = coverage_fitting_image(img, gs)
    @test cov isa Vector{Matrix{Float32}}
    ptexist = CartesianIndex(250, 250)
    K = MMatrix{2, 2, Float64}([80.0 25.0; 10.0 -85.0])
    plot_glyph_given_value!(cov, ptexist, K, gs)
    apply_color_by_coverage!(img, cov, gs)
    mark_at!(img, ptexist; side = 1)
    rc = coarse_radius_for_plotting(gs, K)
    for α = 0:0.005:2π
        # Point on rough outer circle
        Δi = -Int(round(rc * sin(α)))
        Δj = Int(round(rc * cos(α)))
        npt = ptexist + CartesianIndex((Δi, Δj))
        mark_at!(img, npt; side = 1)
        # 
        # From pt to ptexist
        u = ptexist.I .- npt.I
        # The angle for the ray is calculated as if we're in a 'y up' c.s.
        αu = atan(-u[1], u[2]) 
        # αu ≈ α + π
        rptexist = radial_distance_glyph(gs, K, αu + π)
        ptint = ptexist - CartesianIndex(Int.(round.(u .* rptexist ./ norm(u))))
        mark_at!(img, ptint)
    end
    img
    @test is_hash_stored(img, vhash)
end

@testset begin "Tensor glyph crash detection"
    ppts = [CartesianIndex(100, 100), 
            CartesianIndex(100, 110)]
    gs = GSTensor(multip = 0.75, maxg = 100, ming = -100)
    function d!(K, args...)
        K .= [80.0 25.0; 10.0 -85.0]
        K
    end
    d = BidirectionOnGrid(d!, zeros(Float64, 199, 199))
    # Visual check
    img = fill(RGB{N0f8}(0.624,0.667,0.663), 200, 200)
    plot_glyphs!(img, d, ppts, gs)
    # Only the first glyph is accepted (both directions crashing)
    @test placements_and_values(d, gs, ppts) == ([ppts[1]], [[80.0 25.0; 10.0 -85.0]])
    #
    ppts = [CartesianIndex(100, 65), 
            CartesianIndex(90, 155)]
    img = fill(RGB{N0f8}(0.624,0.667,0.663), 200, 200)
    plot_glyphs!(img, d, ppts, gs)
    # Only the first glyph is accepted (primary direction)
    @test length(placements_and_values(d, gs, ppts)[1]) == 1
    #
    img = fill(RGB{N0f8}(0.624,0.667,0.663), 300, 300)
    d = BidirectionOnGrid(d!, zeros(Float64, 299, 299))
    ppts = [CartesianIndex(80, 125), 
            CartesianIndex(217, 160)]
    plot_glyphs!(img, d, ppts, gs)
    # Only the first glyph is accepted (secondary direction crash)
    @test length(placements_and_values(d, gs, ppts)[1]) == 1
    # 
    img = fill(RGB{N0f8}(0.624,0.667,0.663), 300, 300)
    ppts = [CartesianIndex(140, 125), 
            CartesianIndex(140, 185)]
    plot_glyphs!(img, d, ppts, gs)
    # Major axis crossing over minor
    @test length(placements_and_values(d, gs, ppts)[1]) == 1
end


@testset begin "2d vector glyph crash detection, 1 ends within 2"
    vhash = String[]
    COUNT[] = 0
    gs = GSVector(multip = 1.16, maxg = 200)
    # Here, we check if the crash test works with the
    # functor machinery, especially getting coordinate systems
    # right
    function fdir!(v, M)
        α = M[3, 3]
        # Output is interpreted as (x, y) oriented
        v .= [70 * cos(α ), 70 * sin(α)]
        v
    end
    ppts = [CartesianIndex(50, 50), 
            CartesianIndex(50, 150)]
    z = zeros(Float64, 199, 199)
    z[ppts[1]] = -36.6 * π / 180
    z[ppts[2]] = -127 * π / 180
    vog = Vec2OnGrid(fdir!, z)
    @test vog(ppts[1].I...) ≈ [56.19722326337802, -41.73574124759311]
    @test vog(ppts[2].I...) ≈ [-42.127051620643385, -55.90448570331049]
    # Visual check
    img = fill(RGB{N0f8}(0.624,0.667,0.663), 300, 300)
    plot_glyphs!(img, vog, ppts, gs)
    @test is_hash_stored(img, vhash)
    # 
    @test length(placements_and_values(vog, gs, ppts)[1]) == 1
end

@testset begin "2d vector glyph crash detection, 1 crosses 2"
    vhash = String[]
    COUNT[] = 0
    gs = GSVector(multip = 1.16, maxg = 200)
    # Here, we check if the crash test works with the
    # functor machinery, especially getting coordinate systems
    # right
    function fdir!(v, M)
        α = M[3, 3]
        # Output is interpreted as (x, y) oriented
        v .= [70 * cos(α ), 70 * sin(α)]
        v
    end
    ppts = [CartesianIndex(70, 70), 
            CartesianIndex(50, 150)]
    z = zeros(Float64, 199, 199)
    z[ppts[1]] = -36.6 * π / 180
    z[ppts[2]] = -127 * π / 180
    vog = Vec2OnGrid(fdir!, z)
    # Visual check
    img = fill(RGB{N0f8}(0.624,0.667,0.663), 300, 300)
    plot_glyphs!(img, vog, ppts, gs)
    @test is_hash_stored(img, vhash)
    # 
    @test length(placements_and_values(vog, gs, ppts)[1]) == 1
end

@testset begin "Bidirectional vector glyph crash detection"
    vhash = String[]
    COUNT[] = 0
    gs = GSTensor(multip = 1, maxg = 100, ming = -100, direction = 1)
    function fdir!(K, a, b, c, d, M, e, f)
        α = M[3, 3]
        K .= [-85 * sin(α) -85 * cos(α ); 65 * sin(α + 1.5) 65 * cos(α + 1.5)]
        K
    end
    ppts = [CartesianIndex(65, 65), 
            CartesianIndex(100, 150)]
    z = zeros(Float64, 199, 199)
    z[ppts[1]] = -36.6 * π / 180
    z[ppts[2]] = -127 * π / 180
    bdog = BidirectionOnGrid(fdir!, z)
    # 
    img = fill(RGB{N0f8}(0.624,0.667,0.663), 300, 300)
    plot_glyphs!(img, bdog, ppts, gs)
    @test is_hash_stored(img, vhash)
    @test length(placements_and_values(bdog, gs, ppts)[1]) == 1
end

@testset begin "Bidirectional tensor glyph crash detection"
    vhash = String[]
    COUNT[] = 0
    gs = GSTensor(multip = 1, maxg = 100, ming = -100)
    function fdir!(K, a, b, c, d, M, e, f)
        α = M[3, 3]
        K .= [-85 * sin(α) -85 * cos(α ); 65 * sin(α + 1.5) 65 * cos(α + 1.5)]
        K
    end
    ppts = [CartesianIndex(65, 65), 
            CartesianIndex(100, 150)]
    z = zeros(Float64, 199, 199)
    z[ppts[1]] = -36.6 * π / 180
    z[ppts[2]] = -127 * π / 180
    bdog = BidirectionOnGrid(fdir!, z)
    # 
    img = fill(RGB{N0f8}(0.624,0.667,0.663), 300, 300)
    plot_glyphs!(img, bdog, ppts, gs)
    @test is_hash_stored(img, vhash)
    # 
    @test length(placements_and_values(bdog, gs, ppts)[1]) == 1
end
