using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: plot_streamlines!, 𝐊!, 𝐊ᵤ!, Stroke, indices_on_grid
using BitmapMapsExtras: SelectedVec2AtXY
using BitmapMapsExtras: PALETTE_GRGB, MersenneTwister
using BitmapMapsExtras: mark_at!, display_if_vscode
import StatsBase
using StatsBase: Weights, sample

!@isdefined(is_hash_stored) && include("common.jl")

@testset "Streamlines, normalized functor or not" begin
    vhash = String[]
    COUNT[] = 0
    saxy = SelectedVec2AtXY(𝐊!, z_cylinder(π / 6), false, true, normalize =false);
    img = background(saxy)
    pts = indices_on_grid(saxy; Δ = 300)
    mark_at!(img, pts, 5, "in_circle")
    # With non-normalized curvature, 𝐊!, the values du / dt are very small.
    # At time step 1.0, the streamline integration may not have left the initial pixel.
    # That's just a lot of sub-pixel interpolation and extremely in-effective.
    plot_streamlines!(img, saxy, pts; stroke = Stroke(color = PALETTE_GRGB[2]), tstop = 200000)
    # Normalized curvature. The downside is that we can't get an integrated value of curvature, i.e.
    # local inclination, which might be useful for checks but not for finding the streamline.
    saxy = SelectedVec2AtXY(𝐊!, z_cylinder(π / 6), false, true);
    plot_streamlines!(img, saxy, pts; stroke = Stroke(color = PALETTE_GRGB[3], strength = 0.6), dtmax = 1)
    @test is_hash_stored(img, vhash)
end

@testset "Solution keyword example" begin
    vhash = String[]
    COUNT[] = 0
    saxy = SelectedVec2AtXY(𝐊!, z_ellipsoid(), true, true)
    pts = indices_on_grid(saxy)
    img = background(saxy)
    mark_at!(img, pts)
    # This takes too long strides and follow the wrong track
    plot_streamlines!(img, saxy, pts)
    @test is_hash_stored(img, vhash)
    #
    img = background(saxy)
    mark_at!(img, pts)
    # This is persistent
    plot_streamlines!(img, saxy, pts; dtmax = 1.0)
    @test is_hash_stored(img, vhash)
    #
    # Minor direction:
    saxy = SelectedVec2AtXY(𝐊!, z_ellipsoid(), false, true)
    # Also persistent
    plot_streamlines!(img, saxy, pts; dtmax = 1, stroke = Stroke(color = PALETTE_GRGB[3]))
    @test is_hash_stored(img, vhash)
end

@testset "Streamlines starting on grid" begin
    vhash = String[]
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
        z_wavy]
    for fz in vzf
        saxy = SelectedVec2AtXY(𝐊!, fz(), true, true);
        pts = indices_on_grid(saxy)
        img = background(saxy)
        # dtmax = 1.0 is 'always' a good idea.
        plot_streamlines!(img, saxy, pts; dtmax = 1.0)
        display(img)
        #@test is_hash_stored(img, vhash)
    end
end


#@testset "Symmetric directions with different appearance" begin
    vhash = String[]
    COUNT[] = 0
    # Short streamlines from grid points. 
    tstop = 150
    saxy = SelectedVec2AtXY(𝐊!, z_wavy(), false, false)
    pts = indices_on_grid(saxy; Δ = 250, offset = (50, 30))
    pts = [pts[4]]
    img = 0.4N0f8 .* background(saxy) .+ RGB{N0f8}(0.08, 0.08, 0.08)
    plot_streamlines!(img, saxy, pts; dtmax = 1, tstop, dtmin = 0.1)
    mark_at!(img, pts, 5, "in_circle")
    saxy = SelectedVec2AtXY(𝐊!, z_wavy(), false, true)
    stroke = Stroke(color = PALETTE_GRGB[3])
    plot_streamlines!(img, saxy, pts; dtmax = 1, tstop, dtmin = 0.1, stroke)

## TODO: Fix asymmetry. 
    # Short streamlines from grid points. The direction is the second of two.
    @test is_hash_stored(img, vhash)
#end






# Examine point....
using BitmapMapsExtras: paint_convexity!
saxy = SelectedVec2AtXY(𝐊!, z_wavy()[300:300+98, 500:600], false, true);
bdog = saxy.baxy.bid.bdog
img = background(saxy)
pts = CartesianIndex.( [(20, 18), (20, 23)])
pts = [pts[2]]
plot_streamlines!(img, saxy, pts; dtmax = 1, dtmin = 0.1, tstop = 25)
#┌ Debug: In doubt: dotprod = 0.9998989573616911 at (x, y) = [18.442689520680062, 67.59224025282056]
#│ (xprev, yprev) = [18.886201435935202, 68.4884487442526]
x = 18.442689520680062
y = 67.59224025282056
pts = CartesianIndex.( [(Int(floor(saxy.negy(y))), Int(floor(x)))])
mark_at!(img, pts)
pts = CartesianIndex.( [(33, 18)])
mark_at!(img, pts)
x = pts[1][2] + 0.1
y = pts[1][1] + 0.1
saxy(x,y)
mc = saxy.baxy.bid.corners
K1 = mc[1,1]
K2 = mc[2,1]


using BitmapMapsExtras: steepness_radian, z_matrix, norm
steepness_radian(z_matrix(saxy), pts[1]) * 180 / π

# expand dot_product_with_previous
u0 = [12.0, 67.0]
u1 = [11.90773827493078, 66.82227616326072]
Δu = (u1 - u0)
magΔ = norm(Δu)
# ....normalized to unit length
Δu ./= magΔ
# Normalized differential du1 
du = saxy.v
norm(du) < 0.98 && return 1.0
# Dot product    
Δu[1] * du[1] + Δu[2] * du[2]







#plot_glyphs!(img, bdog, pts, GSTensor(;multip = 30))
#
# Confirm conflicting column (major-minor switch)
#
@test bdog(pts[1].I...) == K1
@test bdog(pts[2].I...) == K2
@test is_close_to_perpendicular(K1[:, 1], K2[:, 1])
@test is_close_to_perpendicular(K1[:, 2], K2[:, 2])
# Test the fix by calling the non-discrete abstraction with internal point
di, dj = 0.1, 0.1 
ptint = pts[1].I .+ (di, dj) # This makes pts[1] dominant over pts[2]
x = ptint[2]
y = saxy.negy(ptint[1])
saxy(x,y)
mc = saxy.baxy.bid.corners
@test mc[1,1] == K1
@test mc[2,2][:, 1] == K2[:, 2]
@test mc[2,2][:, 2] == K2[:, 1]
# Now with another dominant_j
di, dj = 0.9, 0.9
ptint = pts[1].I .+ (di, dj) # This makes pts[2] dominant over pts[1]
x = ptint[2]
y = saxy.negy(ptint[1])
saxy(x,y)
mc = saxy.baxy.bid.corners
@test mc[2,2] == K2
@test mc[1,1][:, 1] == K1[:, 2]
@test mc[1,1][:, 2] == K1[:, 1]




























@testset "Many streamlines in both directions from grid points" begin
    vhash = String[]
    COUNT[] = 0
    stroke = Stroke(strength = 0.05)
    tstop = 2000 # default is 1000
    dtmax = 1
    # Major curvature direction. There's a tendency in where they collect,
    # but the tendency is masked by the regular grid
    saxy = SelectedVec2AtXY(𝐊!, z_ridge_peak_valleys(), false, true)
    pts = indices_on_grid(saxy, Δ = 50, offset = 20) # 361
    img = background(saxy; α = 0.4)
    #mark_at!(img, pts)
    plot_streamlines!(img, saxy, pts; dtmax, tstop, stroke)
    saxy = SelectedVec2AtXY(𝐊!, z_ridge_peak_valleys(), false, false)
    plot_streamlines!(img, saxy, pts; dtmax, tstop, stroke)
    @test is_hash_stored(img, vhash)
    #
    # Major curvature direction. There's an other tendency in where they collect,
    # but the tendency is masked by the regular grid
    saxy = SelectedVec2AtXY(𝐊!, z_ridge_peak_valleys(), false, true)
    img = background(saxy; α = 0.4)
    #mark_at!(img, pts)
    plot_streamlines!(img, saxy, pts; dtmax, tstop, stroke)
    saxy = SelectedVec2AtXY(𝐊!, z_ridge_peak_valleys(), false, false)
    plot_streamlines!(img, saxy, pts;  dtmax, tstop, stroke)
    @test is_hash_stored(img, vhash)
end

@testset "Streamlines without grid" begin
    vhash = String[]
    COUNT[] = 0
    stroke = Stroke(strength = 0.05)
    tstop = 2000 # default is 1000
    dtmax = 1
    z = z_ridge_peak_valleys()
    # Uniform random distribution
    n = 361
    n = 10
    pts = sample(MersenneTwister(124), CartesianIndices(z), n)
    # Major curvature. Tends to
    # aggregate on ridges
    img = background(z)
    img .*= 0.2
    img .+= RGB{N0f8}(0.8, 0.8, 0.8)

    plot_streamlines!(img, SelectedVec2AtXY(𝐊!, z, true, true), pts; dtmax, tstop, stroke)
    plot_streamlines!(img, SelectedVec2AtXY(𝐊!, z, true, false), pts; dtmax, tstop, stroke)
    @test is_hash_stored(img, vhash)
    # Minor curvature. These do NOT aggregate on ridges,
    # but seems to identify lines of no curvature, as well as valley bottoms
    img = background(z)
    img .*= 0.2
    img .+= RGB{N0f8}(0.8, 0.8, 0.8)
    plot_streamlines!(img, SelectedVec2AtXY(𝐊!, z, false, true), pts; dtmax, tstop, stroke)
    plot_streamlines!(img, SelectedVec2AtXY(𝐊!, z, false, false), pts; dtmax, tstop, stroke)
    @test is_hash_stored(img, vhash)
end




#= Commented out while BitmapMaps is a troublesome dependency
    z = z_ridge_peak_valleys()

    # Originating at points where z > 0
    seed_dens = clamp.(z, 0.0, 1.0);
    display_if_vscode(background(seed_dens))
    n = 1000
    pts = sample(MersenneTwister(123), CartesianIndices(seed_dens), Weights(vec(seed_dens)), n);
    saxy = SelectedVec2AtXY(𝐊ᵤ!, z, true, true)
    stroke = Stroke(color = PALETTE_GRGB[3], r = 1, strength = 0.3)
    img = background(z)
    mark_at!(img, pts)
    plot_streamlines!(img, saxy, pts; dtmax = 1, stroke)
    saxy = SelectedVec2AtXY(𝐊!, z, true, false)
    plot_streamlines!(img, saxy, pts; dtmax = 1, stroke)
    display_if_vscode(img)
    @test is_hash_stored(img, vhash)
    #
    # Originating at 'sources' (convex terrain)
    vals1 = divergence_of_gradients(-z_ridge_peak_valleys())
    background(vals1)
    vals2 = clamp.(vals1, 0.003, 1.0) .- 0.003
    seed_dens = clamp.(vals2 .* 100, 0, 1.0)
    # Drop the extreme area near the centre.
    seed_dens[400:600, 300:700] .= 0.0
    background(seed_dens)
    n = 1000
    pts = sample(MersenneTwister(123), CartesianIndices(seed_dens), Weights(vec(seed_dens)), n);
    img = background(z)
    mark_at!(img, pts)
    display_if_vscode(img)
    plot_streamlines!(img, saxy, pts; dtmax = 1)
    display_if_vscode(img)
    @test is_hash_stored(img, vhash)
end
=#