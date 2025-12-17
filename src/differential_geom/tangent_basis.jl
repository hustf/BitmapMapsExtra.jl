# This file contains tangent_basis!, tangent_basis  and 
# supporting functions. Mutating versions for tight loops.
# Also 
# - descent!      the normalized surface normal vector projected onto the base plane
# - descent_unit! ...which is the same, but normalized or set to zero.
# All sampled from a 5x5 window.

# Relies on constants KERN1´, KERN2´ and MAG_EPS

"""
    tangent_basis(M)
    ---> P:: 3x3 Matrix

Find the tangent basis for the centre index of M, a 5x5 matrix. M is part of 
a height map or scalar field.

Returns P, the 3×3 matrix whose columns are the tangent basis vectors expressed in the standard basis,
and 
 - where 𝐞₁ lies in the positive standard-basis xy-plane
 - where 𝐞₂ lies on the positive y side
 - where 𝐞₃ lies on the positive z side


P = [𝐞₁ 𝐞₂ 𝐞₃] = 
    [a  d  g;
     0  e  h;
     c  f  i]
"""
function tangent_basis(M)
    # Pre-allocate mutable fixed-sized 2d vector
    v = MVector{2,Float64}(Array{Float64, 1}(undef, 2))
    # Pre-allocate mutable statically-sized array
    P = MMatrix{3,3,Float64}(Array{Float64, 2}(undef, (3, 3)))
    # Mutating P
    tangent_basis!(P, v, M)
end

"""
    tangent_basis!(P, v, M)

Mutating, see `tangent_basis`. Both P and v 
are mutating, v is just an intermediatory to
reduce allocations.
"""
function tangent_basis!(P, v, M)
    # Horizontal direction, r || x, z is out-of screen
    tangent_unit_2d_vector!(v, dz_over_dx, M)
    P[1, 1] = v[1] # r
    P[2, 1] = 0.0
    P[3, 1] = v[2] # z
    # y direction, , r || y, z is out-of screen
    tangent_unit_2d_vector!(v, dz_over_dy, M)
    P[1, 2] = 0.0
    P[2, 2] = v[1] # r
    P[3, 2] = v[2] # z
    # The first column, horizontal direction 𝐞₁ is kept, useful as an angular reference
    # since its y-component is zero.
    #
    # The second column is probably not yet quite orthogonal to the first.
    # We add a vector which is parallel with 𝐞₁ to 𝐞₂. Ref. Gram-Schmidt. After this,
    # there's no component of 𝐞₁ in the second column, but the length of the second column
    # will no longer be 1.
    # 
    # What we are doing:
    # P[:, 2] -= (P[:, 2] ⋅ P[:, 1] ) * P[:, 1]
    # The same, but without array allocation:
    dotval = P[:,2] ⋅ P[:,1] 
    @inbounds for i in 1:3
        P[i,2] -= dotval * P[i,1]
    end
    # Fallback: columns are nearly parallel (degenerate surface).
    # Replace with arbitrary orthonormal vector to maintain valid basis.
    # Tangent meaning is lost in this case.
    if norm(P[:,2]) < 1e-8 # Float64: √1.1×10−16
        _safe_second_vector!(P)
    end
    # Another check: the second column ought to still point screen-up
    if P[2, 2] < 0
        # Flip 180° to point screen up
        P[:, 2] .= -P[:, 2] 
    end
    # (Re-)Adjust the second column to unit length
    normalize!(view(P, :, 2))
    # We'll complete the right-handed orthonormal triad now.
    # This is what we're doing:
    #  P[:, 3] .= P[:, 1] × P[:, 2]
    # The same, but without allocation:
    @inbounds begin
        P[1,3] = P[2,1]*P[3,2] - P[3,1]*P[2,2]
        P[2,3] = P[3,1]*P[1,2] - P[1,1]*P[3,2]
        P[3,3] = P[1,1]*P[2,2] - P[2,1]*P[1,2]
    end
    P
end

function _safe_second_vector!(P)
    println(".")
    throw("Enable if this occurs")
    # Choose an axis that is furthest from P[:, 1]
    if abs(P[1,1]) < abs(P[2,1])
        P[:, 2] .= (1.0, 0.0, 0.0)
    else
        P[:, 2] .= (0.0, 1.0, 0.0)
    end
    P[:,2] .-= (P[:,2] ⋅ P[:,1]) * P[:,1]
    normalize!(P[:,2])
    # remove the component along P[:,1] and normalize
    P[:,2] .= normalize(a - (a⋅P[:,1]) * P[:,1])
end

"""
    tangent_unit_2d_vector(f, M)
    --> (Δr, Δz)

The component Δr = 1 for a flat tangent, and decreases towards 0⁺ for 
steep tangents.

# Arguments

- f   Function, either `dz_over_dy` or `dz_over_dx`.
- M a 5x5 matrix. The tangent vector is estimated for the centre point.
"""
function tangent_unit_2d_vector(f, M)
    # Change over one step forward
    z´ = f(M)
    # The triangle [1, z´] has hypotenus > 1
    l = hypot(1, z´)
    # Scale both sides down to get hypotenus of length 1.
    1 / l, z´ / l
end

"""
    tangent_unit_2d_vector!(v, f, M)

Mutating version of `tangent_unit_2d_vector`, see that.
"""
function tangent_unit_2d_vector!(v, f, M)
    # Change over one step forward
    z´ = f(M)
    # The triangle [1, z´] has hypotenus > 1
    l = hypot(1, z´)
    # Scale both sides down to get hypotenus of length 1.
    v .=  (1 / l, z´ / l)
    v
end


"""
    dz_over_dy(M)
    --> ≈ dz / dy

Smoothed derivative for the center pixel in vertical direction.
Positive direction is upwards if M is part of an image, because
the sign of constant KERN1´ is flipped compared to what is standard.
"""
dz_over_dy(M) = KERN1´ ⋅ M

"""
    dz_over_dx(M)
    --> ≈ dz / dx

Smoothed derivative for the center pixel in horizontal
direction (if M is part of an image).
"""
dz_over_dx(M) = KERN2´ ⋅ M


"""
    angle_xy_to_tangent(a, c, d, e, f, g, h, i, vα::Vector{Float64})
    --> Vector{Float64}
    angle_xy_to_tangent(a, c, d, e, f, g, h, i, α)
    --> Float64
"""
function angle_xy_to_tangent(a, c, d, e, f, g, h, i, vα::Vector{Float64})
    map(α -> angle_xy_to_tangent(a, c, d, e, f, g, h, i, α), vα)
end
function angle_xy_to_tangent(a, c, d, e, f, g, h, i, α)
    @assert i !== 0
    # Symbolic abbreviations that the compiler would probably handle on its own....
    # Might cause unnecessary allocation?
    co = cos(α)
    si = sin(α)
    # 𝐩 lies in the xy-plane (z = 0). We lift it along the z-axis by z,
    # so that the result 𝐫 = 𝐩 + [0, 0, z] lies in the tangent plane.
    # Then we express 𝐫 in the tangent basis and extract its direction.
    z = -(g·co + h·si) / i
    # 𝐫 = 𝐩 + z·𝐤 (in standard basis), then expressed in tangent basis:
    u = a·co +        c⋅z   # component along 𝐞₁
    v = d·co + e·si + f·z  # component along 𝐞₂
    # It sits on a ray from origin.
    ϕ = atan(v, u)
end


"""
    angle_tangent_to_xy(a, d, e, g, h, i, vϕ::Vector)
        --> Vector{Float64}
    angle_tangent_to_xy(a, d, e, g, h, i, ϕ)
        --> Float64
"""
angle_tangent_to_xy(a, d, e, g, h, i, vϕ::Vector) = map(ϕ -> angle_tangent_to_xy(a, d, e, g, h, i, ϕ), vϕ)
function angle_tangent_to_xy(a, d, e, g, h, i, ϕ)
    @assert i !== 0
    # Symbolic abbreviations that the compiler would probably handle on its own....
    # Might cause unnecessary allocation?
    co = cos(ϕ)
    si = sin(ϕ)
    # Tangent-space vector [u, v, w]
    u, v, w = co, si, 0.0
    # Standard-basis vector 𝐩 = P ⋅ 𝐮
    x = a·u + d·v + g·w
    y = e·v + h·w
    # z = c·u + f·v + i·w  # unused
    α = atan(y, x)
end

"""
    descent(M)

Alias: 𝐧ₚ

Projection of surface normal `𝐧` into the xy-plane (y is up, which is a compromise convention here).
Input is a 5x5 matrix of numbers representing elevation. 
Calling `𝐧ₚ!` is preferrable for speed.
"""
function descent(M)
    v = Array{Float64, 1}(undef, 2)
    𝐧ₚ!(v, M)
end



"""
    descent!(v, M)

Alias: 𝐧ₚ!. 

Also see `descent_unit!`.

Projection of `𝐧` into the xy-plane (y is up). `𝐧` is the unit normal vector to the centre of elevation surface `z`.

`M` is a 5x5 view of the elevation surface `z`. This function estimates at `M[3, 3]`.

`v` is mutated to the output.

    v[1]: x ("j")-component (columns) of 𝐧ₚ,
    v[2]: y ("-i")-component of (rows) of 𝐧ₚ 

Refers global const KERN1´  and KERN2´.
KERN2´ values ensure that 'y is up'.

# Example

```
julia> 𝐧ₚ!([0.0, 0.0], [2i for i = 1:5, j = 1:5])
2-element Vector{Float64}:
 -8.068210179770734e-17
  0.894427906538234
```
"""
function descent!(v, M)
    dz_x = dz_over_dx(M)
    dz_y = dz_over_dy(M)
    # Note three components before projection
    mag = hypot(1.0, dz_x, dz_y)
    v[1] = -dz_x / mag
    v[2] = -dz_y / mag
    v
end

"""
    descent_unit!(v, M; minnorm = MAG_EPS)

Alias: 𝐧ₚᵤ!.

See `descent!`.

The 2d vector is normalized to length 1, or set to zero when below threshold `minnorm`. 

Normalized projection of 𝐧 into the xy-plane (y is up). `𝐧` is the normal vector to the elevation surface `z`.
"""
function descent_unit!(v, M; minnorm = MAG_EPS)
    dz_x = dz_over_dx(M)
    dz_y = dz_over_dy(M)
    v[1] = -dz_x
    v[2] = -dz_y
    normalize_or_zero!(view(v, :, 1), minnorm)
end
