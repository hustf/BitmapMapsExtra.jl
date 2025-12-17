# Graphical output of convexity categories 
#
# Relies on constant PALETTE_GRGB

"""
    paint_steepness_rank!(img::Matrix{<:RGB}, z, pts;
        steepness_interval_rad = π / 8, palette = PALETTE_GRGB)
    paint_steepness_rank!(img, z; kws...)

Overlays a color-map onto the image img. Default: 22.5° intervals.

See callees `color_point_by_steepness!` and  `steepness`!
"""
function paint_steepness_rank!(img::Matrix{<:RGB}, z, pts;  
    steepness_interval_rad = π / 8, palette = PALETTE_GRGB) 
    #
    buf = zeros(RGB{N0f8}, size(img)...)
    _paint_steepness_rank!(buf, z, pts, steepness_interval_rad, palette)
    # Take hue and chroma from buf
    chromaticity_over!(img, buf)
    img
end 
function paint_steepness_rank!(img, z; kws...)
    paint_steepness_rank!(img, z, CartesianIndices(img); kws...)
end

function _paint_steepness_rank!(buf, z, pts, steepness_interval_rad, palette)
    # Prepare
    R = CartesianIndices(z)
    Ω = CartesianIndices((-2:2, -2:2))
    minj = R[1][2] - Ω[1][2]
    maxi = R[end][1] - Ω[end][2]
    maxj = R[end][2] - Ω[end][2]
    mini = R[1][1] - Ω[1][2]
    Ri = CartesianIndices((mini:maxi, minj:maxj))
    # Color points one at a time
    for pt in filter(pt -> pt ∈ Ri, sort(vec(pts)))
        steep_rad = steepness_radian(z, pt)
        steepness_rank = Int(steep_rad ÷ steepness_interval_rad + 1)
        set_pixel_from_palette!(buf, pt, steepness_rank, palette)
    end
    buf
end


"""
    steepness_radian(M)
    steepness_radian(z, pt)

Absolute slope in radians. Straightforward and invariant with the choice of coordinates.
See `tangent_basis.jl` for related terms.
"""
function steepness_radian(M)
    @assert size(M) == (5,5)
    dz_x = dz_over_dx(M)
    dz_y = dz_over_dy(M)
    atan(hypot(dz_x, dz_y))
end
function steepness_radian(z, pt)
    Ω = CartesianIndices((-2:2, -2:2))
    win = view(z, Ω .+ pt)
    steepness_radian(win)
end

