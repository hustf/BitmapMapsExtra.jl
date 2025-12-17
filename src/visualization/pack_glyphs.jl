# Smartly spaced glyphs
# These eliminate placements until
# no glyphs overlap.

"""
    coarse_radius_for_plotting(gs::GSTensor{D12, 2}, K)
    coarse_radius_for_plotting(gs::GSTensor{D, 1}, K) where D
    coarse_radius_for_plotting(gs::GSVector, v)

# Example
```
julia> coarse_radius_for_plotting(gs::GSTensor{D12, 2}, K)
```
"""
function coarse_radius_for_plotting(gs::GSTensor{D12, 2}, K)
    # Method for both axes-glyph
    # Extract the half-length of principal directions glyphs
    if is_in_limits(gs, K)
        v1 = @view K[:, 1]
        v2 = @view K[:, 2]
        l1 = norm(v1) * gs.multip
        l2 = norm(v2) * gs.multip
        if is_bidirec_vect_positive(v1)
            # Add min radius
            Δl1 = min(0.5, VECTOR_REL_HALFWIDTH * l1)
        else
            # Add max radius
            Δl1 = max(0.5, VECTOR_REL_HALFWIDTH * l1)
        end
        if is_bidirec_vect_positive(v2)
            # Add min radius
            Δl2 = min(0.5, VECTOR_REL_HALFWIDTH * l2)
        else
            # Add max radius
            Δl2 = max(0.5, VECTOR_REL_HALFWIDTH * l2)
        end
        return max(l1 + Δl1, l2 + Δl2)
    else
        # We return 0, because we do not want
        # this potential glyph to be one of the selected ones.
        return 0.0
    end
end
function coarse_radius_for_plotting(gs::GSTensor{D, 1}, K) where D
    # Method for single-axis glyph
    # Extract the half-length of major and minor principal direction
    # glyphs
    if is_in_limits(gs, K)
        v1 = @view K[:, D]
        l1 = norm(v1) * gs.multip
        if is_bidirec_vect_positive(v1)
            # Add min radius
            Δl = min(0.5, VECTOR_REL_HALFWIDTH * l1)
        else
            # Add max radius
            Δl = max(0.5, VECTOR_REL_HALFWIDTH * l1)
        end
        return l1 + Δl
    else
        # We return 0, because we do not want
        # this potential glyph to be one of the selected ones.
        return 0.0
    end
end
function coarse_radius_for_plotting(gs::GSVector, v)
    if is_in_limits(gs, v)
        norm(v) * gs.multip
    else
        # We return 0, because we do not want
        # this potential glyph to be one of the selected ones.
        return 0.0
    end
end

"""
    mirror_to_right_half(θ)
"""
function mirror_to_right_half(θ)
    #  into [-π, π]
    θ1 = θ > π ?  θ - 2π : θ   
    sign(θ1) * min(abs(θ1), π - abs(θ1))
end

"""
    radial_distance_bidirectional_glyph(::Type{GSVector}, vscaled, α)

See `radial_distance_glyph`.
"""
function radial_distance_bidirectional_glyph(::Type{GSVector}, vscaled, α)
    l = norm(vscaled)
    # Max, min radius
    if is_bidirec_vect_positive(vscaled)
        rA = Float32(max(0.5, VECTOR_REL_HALFWIDTH * l))
        rB = min(0.5f0, rA)
    else
        rB = Float32(max(0.5, VECTOR_REL_HALFWIDTH * l))
        rA = min(0.5f0, rB)
    end
    # The glyph points along vscaled (given in x-y-up directions)
    # Angle of 'pointing axis' rel to x (in x-y-up c.s) 
    β = atan(vscaled[2], vscaled[1])
    # Angle center to tip to side 
    δ = atan(rA - rB, l) 
    # Angle from the glyph direction to α, the line to the other placement.
    # We add 2π where necessary to get positive angles.
    # Since this is a bidirectional, symmetric around γ = π / 2,
    # we're mirroring angles in the 2nd and 3d quadrant
    γ = mod2pi(mirror_to_right_half(α - β))
    # Angle from tip center to tip tangent point
    η = atan(rB, l) 
    @assert rA > 0
    @assert rB > 0
    if 0 < γ <= η
        dist = l + rB
    elseif η < γ <= π / 2 - δ
        # Angle between ray and point on circle rA where the side of the arrow is tangent
        θ = π / 2 - δ - γ
        dist = rA / cos(θ)
    elseif π / 2 - δ < γ <= 3π / 2 + δ
        dist = rA
    elseif 3π / 2 + δ < γ <= 2π - η
        # Angle between ray and point on circle rA where the side of the arrow is tangent
        θ = 3π / 2 + δ - γ
        dist = rA / cos(θ)
    else 
        dist = l + rB
    end
    @assert min(rA, rB) <=  dist
    dist
end
"""
    radial_distance_glyph(gs::GSVector, v, α)
    radial_distance_glyph(gs::GSTensor{D12, 2}, K, α)
    radial_distance_glyph(gs::GSTensor{D, 1}, K, α) where D
    radial_distance_glyph(::Type{GSVector}, vscaled, α)

Distance from a glyph's anchor point along a ray at angle `α` to eclipsing glyph boundary. 
The shape of the glyph is strictly convex.

`α` is the counter-clockwise angle from x in 'x-y up' coordinates between screen horizontal (j or x-axis) and a ray to a point on the circumference.
`v` is the 2-d vector.
`vscaled` is `v * gs.multip`
`K` is a TENSORMAP, i.e. two bidirectional and signed vectors.

Here, we assume that glyph values are within limits.
"""
function radial_distance_glyph(gs::GSVector, v, α)
    radial_distance_glyph(GSVector, gs.multip .* v, α)
end
function radial_distance_glyph(gs::GSTensor{D12, 2}, K, α)
    # Scale the two bidirectional vector glyph values
    vs1 = gs.multip .* view(K, :, 1)
    vs2 = gs.multip .* view(K, :, 2)
    lv1 = radial_distance_bidirectional_glyph(GSVector, vs1, α)
    lv2 = radial_distance_bidirectional_glyph(GSVector, vs2, α)
    # Return the longest radial distance
    return max(lv1, lv2)
end
function radial_distance_glyph(gs::GSTensor{D, 1}, K, α) where D
    # Scale the single bidirectional vector's value
    vscaled = gs.multip .* view(K, :, D)
    return radial_distance_bidirectional_glyph(GSVector, vscaled, α)
end
function radial_distance_glyph(::Type{GSVector}, vscaled, α)
    l = norm(vscaled)
    # Max, min radius
    rA = Float32(max(0.5, VECTOR_REL_HALFWIDTH * l))
    rB = min(0.5f0, rA)
    # The glyph points along vscaled (given in x-y-up directions)
    # Angle of 'pointing axis' rel to x (in x-y-up c.s) 
    β = atan(vscaled[2], vscaled[1]) 
    # Angle center to tip to side 
    δ = atan(rA - rB, l) 
    # Angle from the glyph direction to α, the line to the other placement.
    # We add 2π where necessary to get positive angles.
    γ = mod(α - β, 2π)
    # Angle from tip center to tip tangent point
    η = atan(rB, l) 
    @assert rA > 0
    @assert rB > 0
    if 0 < γ <= η
        dist = l + rB
    elseif η < γ <= π / 2 - δ
        # Angle between ray and point on circle rA where the side of the arrow is tangent
        θ = π / 2 - δ - γ
        dist = rA / cos(θ)
    elseif π / 2 - δ < γ <= 3π / 2 + δ
        dist = rA
    elseif 3π / 2 + δ < γ <= 2π - η
        # Angle between ray and point on circle rA where the side of the arrow is tangent
        θ = 3π / 2 + δ - γ
        dist = rA / cos(θ)
    else 
        dist = l + rB
    end
    @assert min(rA, rB) <=  dist
    dist
end

"""
    too_close(gs::AbstractGlyphSpec, pt1, coarse_r1, v1, pt2, coarse_r2, v2)
    too_close(gs::GSTangentBasis, pt1, pt2)

..If the potential glyph at pt1 would intersect with the
already selected glyph at pt2.

A coarse filtering should be done first. Here, we assume
that glyph values are within limits. Input should be pre-filtered.
"""
function too_close(gs::AbstractGlyphSpec, pt1, coarse_r1, v1, pt2, coarse_r2, v2)
    # Vector from pt1 to pt2
    u = pt2.I .- pt1.I
    # Length of u
    lu = norm(u)
    # Exclusive criterion - fast return
    nooverlap = coarse_r1 + coarse_r2 <= lu  
    nooverlap && return false
    # Vector u's angle to horizontal in x-y-up directions
    αu = atan(-u[1], u[2]) 
    #
    overlap_direct_path(gs, v1, v2, αu, lu) && return true
    #
    # Main axis of glyph 1. v1 is given in (x,y) frame
    αm1 = atan(v1[2], v1[1])
    # Main axis of glyph 2. v2 is given in (x,y) frame
    αm2 = atan(v2[2], v2[1])
    #
    # All other checks works the same way: Find a straight sections
    # path pt1 - pt - pt2. 
    # between anchor points which is entirely within the two glyphs.
    overlap_indirect(gs, pt1, v1, pt2, v2, αu, αm1, αm2)
end
function too_close(gs::GSTangentBasis, pt1, pt2)
    # Vector from pt1 to pt2
    u = pt2.I .- pt1.I
    # Length of u
    lu = norm(u)
    # Criterion
    2√2 * gs.halfsize > lu
end


"""
    overlap_direct_path(gs, v1, v2, αu, lu)
"""
function overlap_direct_path(gs, v1, v2, αu, lu)
    # Intersection dist with pt1's boundary along u
    r1 = radial_distance_glyph(gs, v1, αu)
    # Intersection with glyph 2's boundary
    # A ray from pt2 to pt1 intersects glyph 2's extents at distance
    r2 = radial_distance_glyph(gs, v2, αu + π)
    # The glyphs overlap when
    r1 + r2 >  lu
end

"""
    overlap_indirect(gs::AbstractGlyphSpec, pt1, v1, pt2, v2, αu, αm1, αm2)
"""
function overlap_indirect(gs::AbstractGlyphSpec, pt1, v1, pt2, v2, αu, αm1, αm2)
    axis_intersect_within(gs, pt1, v1, pt2, v2, αm1, αm2) && return true
    exit_of_first_is_in_second(gs, pt1, v1, pt2, v2, αm1, 0) && return true
    exit_of_first_is_in_second(gs, pt2, v2, pt1, v1, αm2, 0) && return true
    # At perp angles to pt1-pt2 anchor points
    exit_of_first_is_in_second(gs, pt1, v1, pt2, v2, αu, π / 2) && return true
    exit_of_first_is_in_second(gs, pt1, v1, pt2, v2, αu, -π / 2) && return true
    exit_of_first_is_in_second(gs, pt2, v2, pt1, v1, αu, π  + π / 2) && return true
    exit_of_first_is_in_second(gs, pt2, v2, pt1, v1, αu, π  - π / 2) && return true
    #
    overlap_indirect_additional_paths(gs::AbstractGlyphSpec, pt1, v1, pt2, v2, αm1, αm2) && return true
    false
end


"""
overlap_indirect_additional_paths(gs, pt1, v1, pt2, v2, αm1, αm2)
overlap_indirect_additional_paths(gs::GSTensor{<:Any, 1}, pt1, v1, pt2, v2, αm1, αm2)
overlap_indirect_additional_paths(gs::GSTensor{D12, 2}, pt1, K1, pt2, K2, αp1, αp2)
"""
overlap_indirect_additional_paths(gs, pt1, v1, pt2, v2, αm1, αm2) = false
function overlap_indirect_additional_paths(gs::GSTensor{<:Any, 1}, pt1, v1, pt2, v2, αm1, αm2)
    exit_of_first_is_in_second(gs, pt1, v1, pt2, v2, αm1, π) && return true
    exit_of_first_is_in_second(gs, pt2, v2, pt1, v1, αm2, π) && return true
    false
end
function overlap_indirect_additional_paths(gs::GSTensor{D12, 2}, pt1, K1, pt2, K2, αp1, αp2)
    exit_of_first_is_in_second(gs, pt1, K1, pt2, K2, αp1, π) && return true
    exit_of_first_is_in_second(gs, pt2, K2, pt1, K1, αp2, π) && return true
    # Prepare for the other checks. `v` is major and `w` is minor direction
    w1 = K1[:, 2] # Glyph 1 in direction αs1
    w2 = K2[:, 2] # Glyph 2 in direction αs2
    # Secondary axes of glyphs. w is given in (x,y) frame
    αs1 = atan(w1[2], w1[1])
    αs2 = atan(w2[2], w2[1])
    exit_of_first_is_in_second(gs, pt1, K1, pt2, K2, αs1, 0) && return true
    exit_of_first_is_in_second(gs, pt2, K2, pt1, K1, αs2, 0) && return true
    exit_of_first_is_in_second(gs, pt1, K1, pt2, K2, αs1, π) && return true
    exit_of_first_is_in_second(gs, pt2, K2, pt1, K1, αs2, π) && return true
    false
end

"""
    axis_intersect_within(gs, pt1, v1, pt2, v2, αp1, αp2)
"""
function axis_intersect_within(gs, pt1, v1, pt2, v2, αp1, αp2)
    if has_two_axes(gs)
        return dual_axes_intersect_within(gs, pt1, v1, pt2, v2, αp1, αp2)
    else
        if !has_external_intersection_point(pt1, pt2, αp1, αp2)
            return false
        end
        i3f, j3f = intersection_floating_point(pt1, pt2, αp1, αp2)
        l1, r1, l2, r2 = tempus_fugit(gs, pt1, v1, pt2, v2, i3f, j3f)
        # Criterion: The path pt1 -> pt2 -> pt3 (l1 + l2) stays within a glyph (r1 + r2) 
        return l1 <= r1 && l2 <= r2
    end
end

"""
    tempus_fugit(gs, pt1, K1, pt2, K2, i3f, j3f)

Length 1-3 
Radial distance from 1 to exit of glyph 1 along 1-3
Length 2-3 
Radial distance from 2 to exit of glyph 2 along 2-3
"""
function tempus_fugit(gs, pt1, K1, pt2, K2, i3f, j3f)
    i1, j1 = pt1.I
    i2, j2 = pt2.I
    # Distances from anchor points to intersection
    l1 = hypot(j3f - j1, i3f - i1)
    l2 = hypot(j3f - j2, i3f - i2)
    # The intersection angle 3 may well be 180° to α
    # Signed angles from anchor to intersection (x-y up) 
    # Note we could perhaps find the angle to intersection without slow trigonometry.
    α1s = atan(-(i3f - i1), j3f - j1)
    α2s = atan(-(i3f - i2), j3f - j2)
    # Distance from pt1's anchor point to exit (the intersection with pt1's glyph boundary)
    r1 = radial_distance_glyph(gs, K1, α1s)
    r2 = radial_distance_glyph(gs, K2, α2s)
    return l1, r1, l2, r2
end

"""
    has_two_axes(gs::AbstractGlyphSpec) = false
    has_two_axes(gs::GSTensor{<:Any, 2}) = true
"""
has_two_axes(gs::AbstractGlyphSpec) = false
has_two_axes(gs::GSTensor{D12, 2}) = true

"""
    has_external_intersection_point(pt1, pt2, α1, α2)
"""
function has_external_intersection_point(pt1, pt2, α1, α2)
    # Image coordinates
    i1, j1 = pt1.I
    i2, j2 = pt2.I
    # Paralell glyphs don't intersect
    α1 ≈ α2 && return false
    # This is an internal intersection
    pt1 == pt2 && return false 
    #
    true
end
function intersection_floating_point(pt1, pt2, α1, α2)
    # Image coordinates
    i1, j1 = pt1.I
    i2, j2 = pt2.I
    # The intersection (floating) point is at 3
    #=
    i3f = i1 - (j3f - j1) * tan(α1)
    i3f = i2 - (j3f - j2) * tan(α2)
    =>
    i1 - (j3f - j1) * tan(α1) = i2 - (j3f - j2) * tan(α2)
    =>
    (-j3f + j1) * tan(α1) = i2 - i1 - (j3f - j2) * tan(α2)
    =>
    -j3f * tan(α1) + j3f * tan(α2) = i2 - i1 + j2 * tan(α2) - j1 * tan(α1)
    =>
    j3f * (tan(α1) - tan(α2)) =  i1 - i2 + j1 * tan(α1) - j2 * tan(α2)
    =>
    j3f  =  (i1 - i2 + j1 * tan(α1) - j2 * tan(α2)) / (tan(α1) - tan(α2))
    =#
    # Input was checked by `has_external_intersection_point`
    j3f  = (i1 - i2 + j1 * tan(α1) - j2 * tan(α2)) / (tan(α1) - tan(α2))
    i3f = i1 - (j3f - j1) * tan(α1)
    i3f, j3f
end 
function dual_axes_intersect_within(gs::GSTensor{D12, 2}, pt1, K1::TENSORMAP, pt2, K2::TENSORMAP, αp1, αp2)
    # Major - major
    if has_external_intersection_point(pt1, pt2, αp1, αp2)
        i3f, j3f = intersection_floating_point(pt1, pt2, αp1, αp2)
        l1, r1, l2, r2 = tempus_fugit(gs, pt1, K1, pt2, K2, i3f, j3f)
        # Criterion: The path pt1 -> pt2 -> pt3 (l1 + l2) stays within a glyph (r1 + r2) 
        if l1 <= r1 && l2 <= r2
            return true
        end
    end
    # Prepare for the other checks. `v` is major and `w` is secondary direction
    w1 = K1[:, 2] # Glyph 1 in direction αs1
    w2 = K2[:, 2] # Glyph 2 in direction αs2
    # Secondary axes of glyphs. w is given in (x,y) frame
    αs1 = atan(w1[2], w1[1])
    αs2 = atan(w2[2], w2[1])
    # Major of 1 - secondary of 2
    if has_external_intersection_point(pt1, pt2, αp1, αs2)
        i3f, j3f = intersection_floating_point(pt1, pt2, αp1, αs2)
        l1, r1, l2, r2 = tempus_fugit(gs, pt1, K1, pt2, K2, i3f, j3f)
        # Criterion: The path pt1 -> pt2 -> pt3 (l1 + l2) stays within a glyph (r1 + r2) 
        if l1 <= r1 && l2 <= r2
            return true
        end
    end
    # Secondary of 1 - primary of 2
    if has_external_intersection_point(pt1, pt2, αs1, αp2)
        i3f, j3f = intersection_floating_point(pt1, pt2, αs1, αp2)
        l1, r1, l2, r2 = tempus_fugit(gs, pt1, K1, pt2, K2, i3f, j3f)
        # Criterion: The path pt1 -> pt2 -> pt3 (l1 + l2) stays within a glyph (r1 + r2) 
        if l1 <= r1 && l2 <= r2
            return true
        end
    end
    # Secondary of 1 - secondary of 2
    if has_external_intersection_point(pt1, pt2, αs1, αs2)
        i3f, j3f = intersection_floating_point(pt1, pt2, αs1, αs2)
        l1, r1, l2, r2 = tempus_fugit(gs, pt1, K1, pt2, K2, i3f, j3f)
        # Criterion: The path pt1 -> pt2 -> pt3 (l1 + l2) stays within a glyph (r1 + r2) 
        if l1 <= r1 && l2 <= r2
            return true
        end
    end
    false
end



"""
     exit_of_first_is_in_second(gs::AbstractGlyphSpec, pt1, v1, pt2, v2, αu, Δαu)

Less simple overlap calculation detection.
Checks if a path from `pt1` to `pt2` is within glyphs 1 or 2.
The path starts a straight line from `pt1` to `exit 1` along direction `αu + Δαu`.
It ends with a straight section from exit to `pt2`.
"""
function exit_of_first_is_in_second(gs::AbstractGlyphSpec, pt1, v1, pt2, v2, αu, Δαu)
    #
    # Angle (in x-y-up directions) of ray from pt1
    α = αu + Δαu
    # Distance from pt1's anchor point to exit (the intersection with pt1's glyph boundary)
    r1 = radial_distance_glyph(gs, v1, α)
    # The exit of the ray from pt1's glyph 
    # Floating point, not CartesianIndex
    i_f = pt1.I[1] - r1 * sin(α)
    j_f = pt1.I[2] + r1 * cos(α)
    # q: Vector from the ray's exit of glyph 1,
    # pointing to pt2's anchor point.
    qi = pt2.I[1] - i_f
    qj = pt2.I[2] - j_f
    # Angle of vector q (in x-y up) c.s.
    β = atan(-qi, qj)
    # Length of vector q
    lq = hypot(qi, qj)
    # The length of a ray from glyph 2's anchor point along q 
    # which is inside of glyph 2.
    r2 = radial_distance_glyph(gs, v2, β + π)
    # Criterion for overlap 
    r2 >  lq
end



"""
    placements_and_values(fij::AbstractIJFunctor, gs::AbstractGlyphSpec, ppts)
    --> (Vector{CartesianIndex{2}}, Vector{typeof(fij(first(ppts).I...)))}

Given a functor `fij` and a glyph spec `gs`, select a subset of `ppts`. These are 
placements for non-overlapping glyphs.
"""
function placements_and_values(fij::AbstractIJFunctor, gs::AbstractGlyphSpec, ppts)
    passed_placements = CartesianIndex{2}[]
    valtyp = typeof(fij(first(ppts).I...))
    # Could be a vector of vectors, or a vector of TENSORMAP
    passed_values = Vector{valtyp}()
    passed_radii = Float64[]
    se = Set(ppts)
    while !isempty(se)
        # Set current point, drop it from the set.
        pt = pop!(se)
        value = fij(pt.I...)
        r = coarse_radius_for_plotting(gs, value)
        if r > MAG_EPS
            if isempty(passed_placements)
                # Only this point to consider: There must be room.
                push!(passed_placements, pt)
                push!(passed_values, copy(value))
                push!(passed_radii, r)
            else
                any_too_close = any(zip(passed_placements, passed_radii, passed_values)) do (ppt, pr, pv)
                    too_close(gs, pt, r, value, ppt, pr, pv)
                end
                # None too close => add this glyph as passed
                if ! any_too_close
                    push!(passed_placements, pt)
                    push!(passed_values, copy(value))
                    push!(passed_radii, r)
                end
            end
        end
    end
    passed_placements, passed_values
end


function placements_only(gs::GSTangentBasis, ppts)
    passed_placements = CartesianIndex{2}[]
    # Could be a vector of vectors, or a vector of TENSORMAP
    se = Set(ppts)
    while !isempty(se)
        # Set current point, drop it from the set.
        pt = pop!(se)
        if isempty(passed_placements)
            # Only this point to consider: There must be room.
            push!(passed_placements, pt)
        else
            any_too_close = any(passed_placements) do ppt
                too_close(gs, pt, ppt)
            end
            # None too close => add this glyph as passed
            if ! any_too_close
                push!(passed_placements, pt)
            end
        end
    end
    passed_placements
end




"""
    pack_glyphs!(img, fij::AbstractIJFunctor, gs::AbstractGlyphSpec; scatterdist = 3.0, seed = MersenneTwister(123))
    pack_glyphs!(img, z::Matrix{<:AbstractFloat}, gs::GSTangentBasis; scatterdist = 3.0, seed = MersenneTwister(123))
    pack_glyphs!(img, z::Matrix{<:AbstractFloat}, gs::AbstractGlyphSpec; scatterdist = 3.0, seed = MersenneTwister(123))
"""
function pack_glyphs!(img, fij::AbstractIJFunctor, gs::AbstractGlyphSpec; scatterdist = 3.0, seed = MersenneTwister(123))
    ppts = indices_scattered(fij; scatterdist, seed)
    filtered_placements, filtered_values = placements_and_values(fij, gs, ppts)
    plot_glyphs_given_values!(img, filtered_placements, filtered_values, gs)
    img
end
function pack_glyphs!(img, z::Matrix{<:AbstractFloat}, gs::GSTangentBasis; scatterdist = 3.0, seed = MersenneTwister(123))
    ppts = indices_scattered(z, gs; scatterdist, seed)
    filtered_placements = placements_only(gs, ppts)
    plot_glyphs!(img, z, filtered_placements, gs)
    img
end
function pack_glyphs!(img, z::Matrix{<:AbstractFloat}, gs::AbstractGlyphSpec; scatterdist = 3.0, seed = MersenneTwister(123))
    pack_glyphs!(img, default_ij_functor(z, gs), gs; scatterdist, seed)
end


