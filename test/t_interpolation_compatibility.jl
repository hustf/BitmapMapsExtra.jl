using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras.TestMatrices: background, z_matrix
using BitmapMapsExtras: SelectedVec2AtXY, 𝐊ᵤ!, 𝐊!, normalize!
using BitmapMapsExtras: is_close_to_opposite, is_close_to_perpendicular
using BitmapMapsExtras: paint_convexity_rank!, RGB, N0f8, mark_at!
using BitmapMapsExtras: GSTensor, plot_glyphs!, MMatrix, BidirectionInDomain, BidirectionOnGrid
using BitmapMapsExtras: paint_steepness_rank!, norm, dot_product, is_close_to_perpendicular

!@isdefined(hash_image) && include("common.jl")

# Identify pairs of neighbouring points with major-minor switch.
# Higher value means more strictly within the criteria
function flip_propensity(bdog, pt)
    pt.I[1] < 3 && return 0.0
    pt.I[2] < 3 && return 0.0
    ptn = pt + CartesianIndex((1,1))
    ptn.I[1] > size(bdog)[1] - 3 && return 0.0
    ptn.I[2] > size(bdog)[2]  - 3 && return 0.0
    v1 = bdog(pt.I...)[:, 1]
    v2 = bdog(ptn.I...)[:, 1]
    #!is_close_to_perpendicular(v1, v2) && return 0.0
    dp1 = abs(dot_product(normalize!(v1), normalize!(v2)))
    dp1 > 0.999 ? 0.0 : 1 - dp1
    #=
    # Also consider the other principal
    v1 = bdog(pt.I...)[:, 2]
    v2 = bdog(ptn.I...)[:, 2]
    !is_close_to_perpendicular(v1, v2) && return 0.0
    dp2 = abs(dot_product(normalize!(v1), normalize!(v2)))
    dp1 + dp2
    =#
end

# Plot glyphs along lines.
# 
"""
    plot_glyphs_along_radial_line!(img, bdog, r0, r1, ψ)

ψ is zero: Down along positive i axis. Rotation around the "image k" axis.
"""
function plot_glyphs_along_radial_line!(img, r0, r1, rngψ, bdog, gs::GSTensor)
    rngr = range(r0, r1, length = 10)
    pts = CartesianIndex{2}[]
    for ψ in rngψ
        for r in rngr
            i = 1 + Int(round(r * cos(ψ)))
            j = 1 + Int(round(r * sin(ψ)))
            if i < 998 && j < 998
                push!(pts, CartesianIndex((i, j)))
            end
        end
    end
    plot_glyphs!(img, bdog, pts, gs)
    img
end

bdog = BidirectionOnGrid(𝐊!, z_sphere_with_bulge())

img = background(bdog)
# Mark high-propensity pixels
for pt in CartesianIndices(size(bdog))
    pr = flip_propensity(bdog, pt)
    if pr > 0.0
        mark_at!(img, pt)
    end
end
# Even though the surface flips at a specific pixel, the sampling
# method is imperfect, so the flip is not that fast. This
# shows the transition distance:
plot_glyphs_along_radial_line!(img, 760, 850, (1:7) * π / 16 , bdog, GSTensor(;multip = 100000, strength = 0.4, ming = -300, maxg = 300 ))





img

mark_at!(img, [CartesianIndex(1,1)], 2 * 760 + 1, "on_circle")
mark_at!(img, [CartesianIndex(1,1)], 2 * 820 + 1, "on_circle")




















pt = CartesianIndex(698, 399)
ptn = pt + CartesianIndex((1,1))
bdog(pt.I...)
bdog(ptn.I...)
# Manually check the neighbour pair and identify the problem.
pt = CartesianIndex(697, 400)
ptn = pt + CartesianIndex((1,1))
bdog(pt.I...)
bdog(ptn.I...)
# This should be a valid pair, too:
plot_glyphs!(img, bdog, [pt], GSTensor(;multip = 52000))
plot_glyphs!(img, bdog, [ptn], GSTensor(;multip = 52000))


    v1 = bdog(pt.I...)[:, 1]
    l1 = norm(v1)
    @show l1
    l1 < 0.0003 
    v1 ./= l1
    v2 = bdog(ptn.I...)[:, 1]
    l2 = norm(v2)
    @show l2
    l2 < 0.0003 
    v2 ./= l1
    @show norm(v1) norm(v2)
    dp1 = abs(dot_product(v1, v2))
    @show dp1
    dp1 > 0.1  # So this is too strict
    # Also consider the other principal
    v1 = bdog(pt.I...)[:, 2]
    l1 = norm(v1)
    @show l1
    l1 < 0.0003 
    v1 ./= l1
    v2 = bdog(ptn.I...)[:, 2]
    l2 = norm(v2)
    @show l2 v2
    l2 < 0.0003 
    v2 ./= l1
    @show norm(v1) norm(v2)
    dp2 = abs(dot_product(v1, v2))
    @show dp2
    dp2 > 0.1  
    0.2 - dp1 - dp2

































#=
#@testset "Interpolation compatibility" begin
    # Overlay curvature lines on convexity rank and elevation

    # Zoom in on region with conflicting bidirections
    z = z_wavy()[840:840+98, 350:450]
    bdog = BidirectionOnGrid(𝐊!, z)
    img = background(bdog)
    pts = CartesianIndex.( [(80, 78), (80, 79)])
    plot_glyphs!(img, bdog, pts, GSTensor(;multip = 30))
    mark_at!(img, pts; side = 7)
    K1 = MMatrix{2, 2, Float64}([6.192627424356922e-5 -0.002699788561524359; 0.009680167472459027 1.0467356583619669e-5])
    K2 =  MMatrix{2, 2, Float64}([-0.00013102941408850384 0.0026862227581094943; 0.009729445992607194 4.2710868775433764e-5])
    #
    # Confirm conflicting minor directions.
    #
    @testset "Confirm direction conflict" begin
        @test bdog(pts[1].I...) == K1
        @test bdog(pts[2].I...) == K2
        @test !is_close_to_opposite(K1[:, 1], K2[:, 1])
        @test  is_close_to_opposite(K1[:, 2], K2[:, 2])
    end
    @testset "Test direction fix" begin
        saxyf = SelectedVec2AtXY(𝐊!, z, false, false)
        # Test the fix by calling the non-discrete abstraction with internal point
        di, dj = 0.0, 0.1 
        # Note that argument order for bdog is "i,j", while saxyf takes (x, y)
        ptint = pts[1].I .+ (di, dj) # This makes pts[1] dominant over pts[2]
        x = ptint[2]
        y = saxyf.negy(ptint[1])
        saxyf(x,y)
        mc = saxyf.baxy.bid.corners
        @test mc[1,1][:, 2] == K1[:, 2]
        @test mc[2,1][:, 2] == K1[:, 2]
        @test mc[1,2][:, 2] == K1[:, 2]
        @test mc[2,2][:, 2] == K1[:, 2]
        @test mc[1,1][:, 1] == K1[:, 1]
        @test mc[2,1][:, 1] !== K1[:, 1]
        @test mc[1,2][:, 1] !== K1[:, 1]
        @test mc[2,2][:, 1] !== K1[:, 1]
        # Now with another dominant_j
        di, dj = 0.0, 0.9
        ptint = pts[1].I .+ (di, dj) # This makes pts[2] dominant over pts[1]
        x = ptint[2]
        y = saxyf.negy(ptint[1])
        saxyf(x,y)
        mc = saxyf.baxy.bid.corners
        @test mc[1,1][:, 2] == K2[:, 2]
        @test mc[2,1][:, 2] == K2[:, 2]
        @test mc[1,2][:, 2] == K2[:, 2]
        @test mc[2,2][:, 2] == K2[:, 2]
        @test mc[1,1][:, 1] == K1[:, 1]
        @test mc[2,1][:, 1] == K1[:, 1]
        @test mc[1,2][:, 1] == K2[:, 1]
        @test mc[2,2][:, 1] == K2[:, 1]
    end
    # Find the most obvious flip point candidates
    function flip_propensity(bdog, pt)
        ptn = pt + CartesianIndex((1,1))
        v1 = bdog(pt.I...)[:, 1]
        l1 = norm(v1)
        @show l1
        l1 < 0.0003 && return 0.2
        v1 ./= l1
        v2 = bdog(ptn.I...)[:, 1]
        l2 = norm(v2)
        @show l2
        l2 < 0.0003 && return 0.2
        v2 ./= l1
        @show norm(v1) norm(v2)
        dp1 = abs(dot_product(v1, v2))
        @show dp1
        #dp1 > 0.1  && return 0.2
        # Also consider the other principal
        v1 = bdog(pt.I...)[:, 2]
        l1 = norm(v1)
        @show l1
        l1 < 0.0003 && return 0.2
        v1 ./= l1
        v2 = bdog(ptn.I...)[:, 2]
        l2 = norm(v2)
        @show l2 v2
        l2 < 0.0003 && return 0.2
        v2 ./= l1
        @show norm(v1) norm(v2)
        dp2 = abs(dot_product(v1, v2))
        @show dp2
        dp2 > 0.1  && return 0.2
        0.2 - dp1 - dp2
    end

principal_curvatures_paraboloid 
z_sphere_with_bulge
z_cos                            
z_cylinder                       
z_cylinder_offset                
z_ellipsoid                      
z_exp3                           
z_paraboloid                     
z_plane                          
z_ridge_peak_valleys             
z_sphere                         
z_sphere_with_bulge

   # Zoom in on region with conflicting column order (major-minor switch)
    saxy = SelectedVec2AtXY(𝐊!, z_sphere_with_bulge(), true, false);
    bdog = saxy.baxy.bid.bdog
    background(bdog)
    img = paint_convexity(bdog)
    pts = CartesianIndices((20:880, 20:880))
    for pt in pts
        pr = flip_propensity(bdog, pt)
        if pr < 0.2
            mark_at!(img, pt)
            print(pt)
        end
    end
    img

    background(mdp)
    minimum(mdp)
    findall(<(0.1), mdp) # 29,37
    pts = CartesianIndex.( [(29, 37), (30, 38)])
    img = background(bdog)
    img = paint_convexity(bdog)
    plot_glyphs!(img, bdog, pts, GSTensor(;multip = 3000))
    K1 = MMatrix{2, 2, Float64}([-0.000805262817756826 -0.005899380252591684; 0.00046577943554938816 -0.010243478026677622])
    K2 = MMatrix{2, 2, Float64}([-0.0005206967775479895 -0.006150383329186725; 0.00027887526545671075 -0.011349416886604008])
    #
    # Confirm conflicting column (major-minor switch)
    #
    @testset "Confirm direction conflict" begin
        @test bdog(pts[1].I...) == K1
        @test bdog(pts[2].I...) == K2
        @test is_close_to_perpendicular(K1[:, 1], K2[:, 1])
        @test is_close_to_perpendicular(K1[:, 2], K2[:, 2])
    end
    @testset "Test major minor switch" begin
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
    end
#end

    # Small region with flat-flat curvature at the top,
    # and a flat-convex curvature at the bottom.
    # Axes aligned with global axes
    z = z_cylinder(0.1)[500:500+98, 75:175]
    for i in 38:-1:2
        z[i - 1, :] .= z[i, :]
    end
    # We're using the non-normalized version here, because
    # relative values are of interest
    saxy = SelectedVec2AtXY(𝐊!, z, true, false);
    bdog = saxy.baxy.bid.bdog
    img = paint_convexity(bdog) #, flatval = 0.0001)
    pts = CartesianIndex.( [(35, 50), (36, 50), (38, 50), (39, 50), (40, 50), (40, 50)])
    img[pts] .= RGB{N0f8}(0, 0, 0); img
    plot_glyphs!(img, bdog, pts, GSTensor(;multip = 10000))
    K1 = MMatrix{2, 2, Float64}([0.0 -2.004079560501466e-5; 0.0 -1.0950221588857005e-14])
    K2 = MMatrix{2, 2, Float64}([0.0 -2.004079560501466e-5; 0.0 -1.0950221588857005e-14])    
    K3 = MMatrix{2, 2, Float64}([4.117929205109241e-5 9.137724698591841e-5; 1.299516517375901e-5 -0.0002895577180704922])
    K4 = MMatrix{2, 2, Float64}([-9.24934562373533e-5 0.00020713532313113286; -1.3776187003279276e-5 -0.0013907081706103417])
    K5 = MMatrix{2, 2, Float64}([-9.14028138289192e-5 0.00020977363451338632; -9.608339921479404e-6 -0.00199554709558102])
    K6 = MMatrix{2, 2, Float64}([-9.14028138289192e-5 0.00020977363451338632; -9.608339921479404e-6 -0.00199554709558102])

    #
    # Confirm zero
    #
    #@testset "Confirm direction conflict" begin
        @test bdog(pts[1].I...) == K1
        @test bdog(pts[2].I...) == K2
        @test bdog(pts[3].I...) == K3
        @test bdog(pts[4].I...) == K4
        @test bdog(pts[5].I...) == K5
        @test bdog(pts[6].I...) == K6

        @test norm(K1[:, 1]) < bdog.minnorm
        @test norm(K2[:, 1]) < bdog.minnorm
        @test norm(K3[:, 1]) > bdog.minnorm
        @test norm(K4[:, 1]) > bdog.minnorm
        @test norm(K5[:, 1]) > bdog.minnorm
        @test norm(K6[:, 1]) > bdog.minnorm

        @test norm(K1[:, 2]) > bdog.minnorm
        @test norm(K2[:, 2]) > bdog.minnorm
        @test norm(K3[:, 2]) > bdog.minnorm
        @test norm(K4[:, 2]) > bdog.minnorm
        @test norm(K5[:, 2]) > bdog.minnorm
        @test norm(K6[:, 2]) > bdog.minnorm


    #end

        # Test the fix by calling the non-discrete abstraction with internal point
        di, dj = 0.1, 0.1
        ptint = pts[2].I .+ (di, dj) # This makes pts[1] dominant over pts[2]
        x = ptint[2]
        y = saxy.negy(ptint[1])
        saxy(x, y)
        mc = saxy.baxy.bid.corners
        @test mc[1,1][:, 1] == K2[:, 1]
        @test mc[2,1][:, 1] == K2[:, 1]
        @test mc[1,2][:, 1] == K2[:, 1]
        @test mc[2,2][:, 1] == K2[:, 1]

    end
=#