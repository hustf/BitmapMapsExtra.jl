# Graphical output of convexity categories 
#
# Relies on constant PALETTE_CONVEXITY

"""
    paint_convexity_rank!(img::Matrix{<:RGB}, z::AbstractMatrix{<:Real}, pts;  
        flatval = 0.0003, palette = PALETTE_CONVEXITY)
    paint_convexity_rank!(img, f::T, pts; kws...) where T<: Union{AbstractXYFunctor,
            AbstractIJFunctor, BidirectionAtXY, BidirectionInDomain,
            Vec2InDomain}
    --> img::Matrix{<:RGB}

Overlays a color-map onto the image img, preserving lightness.
See `convexity_rank`!
"""
function paint_convexity_rank!(img::Matrix{<:RGB}, z::AbstractMatrix{<:Real}, pts;  
    flatval = 0.0003, palette = PALETTE_CONVEXITY) 
    #
    buf = zeros(RGB{N0f8}, size(img)...)
    _paint_convexity_rank!(buf, z, pts, flatval, palette)
    # Take hue and chroma from buf
    chromaticity_over!(img, buf)
    img
end 
function paint_convexity_rank!(img, f::T, pts; kws...) where T<: Union{AbstractXYFunctor,
            AbstractIJFunctor, BidirectionAtXY, BidirectionInDomain,
            Vec2InDomain}
    paint_convexity_rank!(img, z_matrix(f), pts; kws...)
end

function _paint_convexity_rank!(buf, z, pts, flatval, palette)
    # Prepare
    bog = BidirectionOnGrid(𝐊!, z)
    R = CartesianIndices(z)
    Ω = CartesianIndices((-2:2, -2:2))
    minj = R[1][2] - Ω[1][2]
    maxi = R[end][1] - Ω[end][2]
    maxj = R[end][2] - Ω[end][2]
    mini = R[1][1] - Ω[1][2]
    Ri = CartesianIndices((mini:maxi, minj:maxj))
    # Color points one at a time
    for pt in filter(pt -> pt ∈ Ri, sort(vec(pts)))
        # K is a 2x2 matrix with two bidirectional vectors.
        K = bog(pt.I...)
        rank = convexity_rank(K, flatval)
        set_pixel_from_palette!(buf, pt, rank, palette)
    end
    buf
end
