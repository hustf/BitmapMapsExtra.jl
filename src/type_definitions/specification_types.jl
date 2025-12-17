 # This defines types for
 # specifying 
 # - glyph properties.
 # - streamline properties
abstract type AbstractGlyphSpec end 

"""
    GSTensor(; multip = 50, ming = -50, maxg = 50,
                  dashsize = 5.0f0, strength = 0.8f0,
                  colors = (PALETTE_GRGB[3], PALETTE_GRGB[4]),
                  direction = 1:2)

Specifies how tensor map glyphs are presented.
"""
struct GSTensor{D, N} <: AbstractGlyphSpec
    multip::Float64
    ming::Float64
    maxg::Float64
    dashsize::Float32
    strength::Float32 # Graphical spray strength
    colors::NTuple{N, RGB{N0f8}}
    function GSTensor{D,N}(multip, ming, maxg, dashsize, strength, colors) where {D,N}
        @assert D == D1 || D == D2 || D == D12
        @assert N == ((D == D12) ? 2 : 1)
        new{D,N}(multip, ming, maxg, dashsize, strength, colors)
    end
end
# Constructor
function GSTensor(; multip = 50, ming = -50, maxg = 50,
                  dashsize = 5.0f0, strength = 0.8f0,
                  colors = (PALETTE_GRGB[3], PALETTE_GRGB[4]),
                  direction = 1:2)
    D, colrs = _dirs_and_colors(direction, colors...)
    N = length(colrs)  # 1 or 2
    GSTensor{D,N}(float(multip), float(ming), float(maxg),
                  Float32(dashsize), Float32(strength), colrs)
end
_dirs_and_colors(d::Int, c1, c2) = d == 1 ? (D1,  (c1,)) :
                                   d == 2 ? (D2,  (c2,)) :
                                            error("direction must be 1, 2, or 1:2")
function _dirs_and_colors(rd::UnitRange{Int}, c1, c2)
    f, l = first(rd), last(rd)
    if     f == 1 && l == 1
        return (D1,  (c1,))
    elseif f == 2 && l == 2
        return (D2,  (c2,))
    elseif f == 1 && l == 2
        return (D12, (c1, c2))
    else
        error("direction must be 1, 2, or 1:2 (got $rd)")
    end
end


@define_show_with_fieldnames GSTensor
function is_in_limits(gs::GSTensor{D12, 2}, K::TENSORMAP)
    K1 = @view K[:, 1]
    K2 = @view K[:, 2]
    limit1 = is_bidirec_vect_positive(K1) ? gs.maxg : abs(gs.ming)
    limit2 = is_bidirec_vect_positive(K2) ? gs.maxg : abs(gs.ming)
    norm(K1) * gs.multip ≤ limit1 || return false
    norm(K2) * gs.multip ≤ limit2 || return false
    return true
end
function is_in_limits(gs::GSTensor{D1, 1}, K::TENSORMAP)
    K1 = @view K[:, 1]
    limit1 = is_bidirec_vect_positive(K1) ? gs.maxg : abs(gs.ming)
    norm(K1) * gs.multip ≤ limit1 || return false
    return true
end
function is_in_limits(gs::GSTensor{D2, 1}, K::TENSORMAP)
    K2 = @view K[:, 2]
    limit2 = is_bidirec_vect_positive(K2) ? gs.maxg : abs(gs.ming)
    norm(K2) * gs.multip ≤ limit2 || return false
    return true
end



"""
    GSVector(; multip = 50, ming = 0.01, maxg = 50,
        dashsize = 5.0f0, strength = 0.7f0, color = COLOR_CURVGLYPH)

Specifies how vector glyphs are presented.
"""
struct GSVector <: AbstractGlyphSpec
    multip::Float64
    ming::Float64
    maxg::Float64
    dashsize::Float32
    strength::Float32 # Graphical spray strength
    color::RGB{N0f8}
end
function GSVector(; multip = 50, ming = 0.01, maxg = 50,
    dashsize = 5.0f0, strength = 0.7f0, color = COLOR_CURVGLYPH)
    #
    GSVector(float(multip), float(ming), float(maxg), Float32(dashsize), Float32(strength), color)
end
@define_show_with_fieldnames GSVector
is_in_limits(gs::GSVector, v) = gs.ming ≤ gs.multip * norm(v) ≤ gs.maxg


struct GSTangentBasis <: AbstractGlyphSpec
    halfsize::Int64
    colors::NTuple{3, RGB{N0f8}}
end
GSTangentBasis(;halfsize = 30, colors = PALETTE_RGB) = GSTangentBasis(halfsize, colors)
@define_show_with_fieldnames GSTangentBasis


abstract type AbstractCurveSpec end
"""
Stroke(;r = 1.2f0, strength = 0.7f0, color = PALETTE_GRGB[4], tdensity = 0.05)

- `tdensity` - How many points per "time" unit we extract from the differential equation's continuous solution. 
- `r` - Spray radius. This applies after the touched pixels are identified. Around each pixel, apply a 
   sprayed circle. Radius > 3 => the centre pixels are not sprayed (since they're probably covered by neighbouring sprays).
- `strength` The nominal coverage applied in each spray (the coverage tapers off). 
   Coverage accumulates 'infinitely' if a pixel is sprayed 'infinitely' many times. After all 
    spraying along streamlines is finished, coverage will be converted non-linearly to color 
    application, see `apply_color_by_coverage`. Hence, small `strength` makes streamlines less 
    visible, unless where many streamlines overlap.
- `color`
"""
struct Stroke <: AbstractCurveSpec
    r::Float32
    strength::Float32 # Graphical spray strength
    color::RGB{N0f8}
    tdensity::Float64 # Points per "time"
end
function Stroke(;r = 1.2f0, strength = 0.7f0, color = PALETTE_GRGB[4], tdensity = 0.05)
    Stroke(Float32(r), Float32(strength), color, tdensity)
end
@define_show_with_fieldnames Stroke