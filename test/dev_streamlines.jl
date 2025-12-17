using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras.TestMatrices: I0
using BitmapMapsExtras: plot_streamlines!, PALETTE_GRGB
using BitmapMapsExtras: SelectedVec2AtXY, 𝐊!, 𝐧ₚ!, Vec2AtXY
using BitmapMapsExtras: Stroke, plot_streamlines, allocations_curvature
using BitmapMapsExtras: is_close_to_perpendicular
using BitmapMapsExtras: RGB, N0f8, paint_convexity_rank!, indices_on_grid, 𝐊ᵤ!
using BitmapMapsExtras: mark_at!, get_streamlines_xy, DrawAndSpray
using BitmapMapsExtras: get_streamlines_points, is_bidirec_vect_positive
using BitmapMapsExtras: dot_product_with_previous
using BitmapMapsExtras.DrawAndSpray: apply_color_by_coverage!
using BenchmarkTools

!@isdefined(hash_image) && include("common.jl")

#r = TestMatrices.r
pts = [CartesianIndex(300, 450),
      CartesianIndex(300, 475),
      CartesianIndex(300, 500)]
z = z_ellipsoid(; tilt = -π / 5 );
img = background(z)
stroke = Stroke(r = 1f0)
saxy = SelectedVec2AtXY(𝐊!, z, true, true);
vaxy = Vec2AtXY(𝐧ₚ!, z);

plot_streamlines!(img, saxy, pts; stroke, dtmax = 1)
plot_streamlines!(img, vaxy, pts; stroke = Stroke(color = PALETTE_GRGB[2]), dtmax = 1)
saxy = SelectedVec2AtXY(𝐊!, z, true, false);
plot_streamlines!(img, saxy, pts; stroke, dtmax = 1)
saxy = SelectedVec2AtXY(𝐊!, z, false, true);
plot_streamlines!(img, saxy, pts; stroke, dtmax = 1)
saxy = SelectedVec2AtXY(𝐊!, z, false, false);
plot_streamlines!(img, saxy, pts; stroke, dtmax = 1)

# Optimize is_close_to_perpendicular / is_pointing_roughly_perpendicular

Ri, Ω, v, P, K, vα, vκ, vβ, lpc = allocations_curvature(TestMatrices.R)

K .= [1 0; 0.1 1]
# 63.367 ns (2 allocations: 64 bytes)
# 79.835 ns (2 allocations: 64 bytes)
@btime is_close_to_perpendicular(K[:, 1], K[:, 2])
K .= [1 0; 0.1 1] 
# 94.737 ns (2 allocations: 64 bytes)
#112.882 ns (2 allocations: 64 bytes)
#92.767 ns (2 allocations: 64 bytes)
# @btime is_pointing_roughly_perpendicular(saxy, K[:, 1], K[:, 2])


function plot_convexity_and_lines_of_curvature(z, major::Bool; pts = indices_on_grid(z; offset = (20, 20)), odekws...)
    saxyf = SelectedVec2AtXY(𝐊ᵤ!, z, major, true);
    saxyb = SelectedVec2AtXY(𝐊ᵤ!, z, major, false);
    img = paint_convexity(saxyf.baxy.bid.bdog.z; flatval = 0.0003, bckpar1 = 0.5f0, bckpar2 = 0.5f0)    
    stroke = Stroke(; color = PALETTE_GRGB[4], strength = 0.4)
    plot_streamlines!(img, saxyf, pts; stroke, dtmax = 1.0, odekws...)
    #stroke = Stroke(color = PALETTE_GRGB[2], strength = 1.0)
    plot_streamlines!(img, saxyb, pts; stroke, dtmax = 1.0, odekws...)
    mark_at!(img, pts)
end

# Explore long curvature streamlines

function streamline_lengths(saxy; pts = indices_on_grid(saxy.baxy.bid.bdog.z), odekws...)
    sols = get_streamlines_points(saxy, pts, 0.05; odekws..., dtmax = 1.0)
    Int.(round.(length.(sols)))
end


function streamline_lengths(saxyf, saxyb; pts = indices_on_grid(saxyf.baxy.bid.bdog.z), odekws...)
    # Provided normalized and dtmax = 1, we can skip point extraction, which is costly.
    # Fwd
    vlb = streamline_lengths(saxyb; pts, odekws...)
    # Bck
    vlf = streamline_lengths(saxyf; pts, odekws...)
    vlb .+ vlf
end

function streamline_lengths(z::Matrix, major, pts; odekws...)
    saxyf = SelectedVec2AtXY(𝐊ᵤ!, z, major, true)
    saxyb = SelectedVec2AtXY(𝐊ᵤ!, z, major, false)
    streamline_lengths(saxyf, saxyb; pts, odekws...)
end

#pts = indices_on_grid(TestMatrices.R; Δ = 80, offset = (10, 10))
pts = CartesianIndex.([(795, 280),
                       (800, 280),
                       (805, 280),
                       (810, 280),
                       (820, 280),
                       (830, 280),
                       (840, 280),
                       (850, 280),
                       (860, 280),
                       (870, 280),
                       (880, 280)])
plot_convexity_and_lines_of_curvature(z_wavy(), false; pts)#, tstop = 500)




pts = CartesianIndex.([(680, 380),
                       (690, 370),
                       (700, 360),
                       (710, 350),
                       (720, 340),
                       (730, 330),
                       (740, 320),
                       (750, 310),
                       (760, 300),
                       (770, 290)])


plot_convexity_and_lines_of_curvature(z_ellipsoid(), true; pts, tstop = 1300)
plot_convexity_and_lines_of_curvature(z_ellipsoid(), false; pts, tstop = 1300)
pts = CartesianIndex.([( 380, 680),
                       ( 370, 690),
                       ( 360, 700),
                       ( 350, 710),
                       ( 340, 720),
                       ( 330, 730),
                       ( 320, 740),
                       ( 310, 750),
                       ( 300, 760),
                       ( 290, 770)])

plot_convexity_and_lines_of_curvature(z_sphere_with_bulge(), false; pts, tstop = 1300)

# Steepness is not the best criterion.
# We should have some indication of when the curvatures are too equal to determine
# direction in a thrustworthy way, and then terminate.



pts = CartesianIndex.([(875, 501)])
pt = pts[1]
saxyf = SelectedVec2AtXY(𝐊ᵤ!, z_wavy(), false, true);
y = float(saxyf.negy(pt.I[1]))
x = float(pt.I[2])
saxyf(x, y) # right
saxyf(x + 1, y) # right
saxyf(x + 2, y) # left
u0 = [x + 1, y]
u1 = [x + 2, y]
dot_product_with_previous(saxyf, u0, u1)
BitmapMapsExtras.is_close_to_opposite(saxyf, u0, u1)
is_bidirec_vect_positive(saxyf(x + 1, y) )
is_bidirec_vect_positive(saxyf(x + 2, y) )
img = background(saxyf)
plot_streamlines!(img, saxyf, pts; dtmax = 1.0, stroke = Stroke(;r = 3f0))
saxyf(x + 1.2, y) 
saxyf(x + 1, y) # interpolation excemption NOT working properly. Beware incorrect comparison '==' 
saxyf(x + 1.6, y) # interpolation excemption working properly.

@edit saxyf(x + 1.5, y)
saxyf.flip[]
x = 502.5
saxyf.baxy(x, y)*-
baxy = saxyf.baxy
baxy.bid(x, baxy.negy(y)) # interpolation excemption not working
bid = baxy.bid
BitmapMapsExtras.update_corners!(bid, i1, j1, i2, j2)
corners = bid.corners
 x - j1
negy - i1
BitmapMapsExtras.interpolate_and_normalize_directions!(bid.lastvalue, bid.corners, x - j1, negy - i1)

    f11 = corners[1, 1] # (x, y) = (0, 0) 
    f12 = corners[1, 2] # (x, y) = (0, 1)
    f21 = corners[2, 1] # (x, y) = (1, 0)
    f22 = corners[2, 2] # (x, y) = (1, 1)



vl = streamline_lengths(z_ridge_peak_valleys(), true, pts; dtmax = 1, tstop = 15000, dtmin = 0.1)
maximum(vl)
sort(vl)
# Visualize where the long ones are
makeodd(x) = isodd(x) ? x : x + 1
function markit!(cov, pt, l)
    r = Int(round(0.05 * l))
    mark_at!(cov, pt; side = makeodd(r), f_is_filled = DrawAndSpray.func_is_in_square)
end
markthose!(cov, pts, vl) = map( zip(pts, vl)) do (pt, l)
    markit!(cov, pt, l)
end
img = paint_convexityloc(z_ridge_peak_valleys())
cov = falses(size(img)...)
markthose!(cov, pts, vl)
sum(cov) / length(cov)
display_if_vscode(cov)
apply_color_by_coverage!(img, 0.8f0 .* cov, RGB{N0f8}(0.9, 0.9, 0.9))

# Let's plot the longest only.
inds = findall(i -> vl[i] > 1836,  1:length(vl))
vl[inds]
img1 = plot_convexity_and_lines_of_curvature(z_ridge_peak_valleys(), true; pts = [pts[inds[1]]], tstop = 55000, dtmin = 0.1)
img2 = plot_convexity_and_lines_of_curvature(z_ridge_peak_valleys(), true; pts = [pts[inds[2]]], tstop = 55000, dtmin = 0.1)
img3 = plot_convexity_and_lines_of_curvature(z_ridge_peak_valleys(), true; pts = [pts[inds[3]]], tstop = 55000, dtmin = 0.1)




# ahh... these may seem to loop back on itself. But it's really a bug in the point extraction....
# It seems to have issues with zigzagging solutions or something.
saxyf = SelectedVec2AtXY(𝐊ᵤ!, z_ridge_peak_valleys(), true, true)
solpts = get_streamlines_points(saxyf, [pts[inds[2]]], 0.05; tstop = 555000, dtmax = 1.0)

length(unique(solpts[1])) < length(solpts[1])
using UnicodePlots

lineplot(getindex.(Tuple.(solpts[1]), 1))
lineplot(getindex.(Tuple.(solpts[1]), 2))

# Dive into get_streamlines_points.
sols = get_streamlines_xy(saxyf, [pts[inds[2]]]; tstop = 555000, dtmax = 1.0)
negy = saxyf.negy
sol_density = 0.05
sol = sols[1]
length(sol.u)
# Dive into extract_discrete_points_on_streamline
ptss = CartesianIndex{2}[]
oldi, oldj = 0, 0
oldoldi, oldoldj = 0, 0
told = 0.0

step = sol_density * sign( last(sol.t) - first(sol.t)   )
# Just the last seconds
t0 = last(sol.t) - 20
t1 = last(sol.t)
trng = range(t0, t1; step)
length(trng)
for t in trng
    x, y = sol(t) # This is a fast, interpolated lookup 
    ny = negy(y)
    i = Int(round(ny))
    j = Int(round(x))
    if i !== oldi || j !== oldj
#        if i !== oldoldi || j !== oldoldi
            # We didn't just visit this pixel before.
            p = CartesianIndex(i, j)
            if p ∈ ptss
                printstyled((i,j), t - told, "\n", color = :yellow)
            else
                printstyled((i,j), t - told, "\n", color = :blue)
            end
            push!(ptss, p)
            oldoldi, oldoldj = oldi, oldj
            oldi, oldj = i, j
            told = t
#        end
    end
end
length(ptss) # 13
length(unique(ptss))
Tuple.(ptss)
# We would need to check a lot of previous time steps to avoid this duplication of 
# points. And adding to a set is hardly a good solution, either (what if we HAVE long loops?).
# That would only mask what is occuring. Setting dtmin seems to be a better solution here.
img
mark_at!(img, CartesianIndex((87,284)), side = 11, f_is_filled = DrawAndSpray.func_is_in_square)
# It shouldn't be a surprise that this oscillation occurs on a ridge.
# Is setting dtmin a good solution everywhere?











# Now major principal
vl = streamline_lengths(z_ridge_peak_valleys(), true, pts; dtmax = 1)
img = paint_convexityloc(z_ridge_peak_valleys())
cov = falses(size(img)...);
markthose!(cov, pts, vl);
apply_color_by_coverage!(img, 0.8f0 .* cov, RGB{N0f8}(0.9, 0.9, 0.9))
img1 = plot_convexity_and_lines_of_curvature(z_ridge_peak_valleys(), true; pts)

img


# Seems we did a bad interpretation, or we discovered a problem.
# Those lines don't seem long at all.

function investigate_sols(z::Matrix, major, pts; odekws...)
    saxyf = SelectedVec2AtXY(𝐊ᵤ!, z, major, true)
    saxyb = SelectedVec2AtXY(𝐊ᵤ!, z, major, false)
    sols1 = get_streamlines_xy(saxyf, pts; odekws...)
    sols2 = get_streamlines_xy(saxyb, pts; odekws...)
    sols1, sols2
end
sols1, sols2 = investigate_sols(z_ridge_peak_valleys(), true, pts[inds])

# check
length(sols1.u[1]) + length(sols2.u[1])
length(sols1.u[2]) + length(sols2.u[2])

# So what happens?

sols1.u[2] # Many of these are close. dtmax?