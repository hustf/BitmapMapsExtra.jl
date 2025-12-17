"""
This module is included in module BitmapMapsExtras.

# Example

```
julia> using BitmapMapsExtras

julia> using BitmapMapsExtras.TestMatrices

julia> using BitmapMapsExtras.TestMatrices: I0

julia> varinfo(TestMatrices)
  name                                  size summary
  ––––––––––––––––––––––––––––––– –––––––––– –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
  TestMatrices                    33.809 KiB Module
  background                         0 bytes background (generic function with 2 methods)
  principal_curvatures_paraboloid    0 bytes principal_curvatures_paraboloid (generic function with 2 methods)
  z_sphere_with_bulge                0 bytes z_sphere_with_bulge (generic function with 1 method)
  z_cos                              0 bytes z_cos (generic function with 1 method)
  z_cylinder                         0 bytes z_cylinder (generic function with 1 method)
  z_cylinder_offset                  0 bytes z_cylinder_offset (generic function with 1 method)
  z_ellipsoid                        0 bytes z_ellipsoid (generic function with 1 method)
  z_exp3                             0 bytes z_exp3 (generic function with 1 method)
  z_paraboloid                       0 bytes z_paraboloid (generic function with 1 method)
  z_plane                            0 bytes z_plane (generic function with 1 method)
  z_ridge_peak_valleys               0 bytes z_ridge_peak_valleys (generic function with 1 method)
  z_sphere                           0 bytes z_sphere (generic function with 1 method)
  z_wavy                             0 bytes z_wavy (generic function with 1 method)
```

"""
module TestMatrices
import ..SelectedVec2AtXY, ..Vec2AtXY, ..z_matrix, ..PALETTE_BACKGROUND
import ..AbstractIJFunctor, ..AbstractXYFunctor
import ..BidirectionAtXY, ..BidirectionInDomain, ..Vec2InDomain
using ..NonlinearSolve: IntervalNonlinearProblem
using ..SciMLBase: successful_retcode, solve
import ImageCore
using ImageCore: scaleminmax, RGB, N0f8
export z_cylinder, z_cylinder_offset, z_sphere, z_ellipsoid, 
    z_paraboloid, z_cos, z_exp3, z_plane, z_ridge_peak_valleys,
    z_wavy, z_sphere_with_bulge
export principal_curvatures_saddle, principal_curvatures_paraboloid
export background
##############################################
# Define (parameters for) test array functions
##############################################
const w_h = 999
const r = w_h ÷ 2
const R = CartesianIndices((w_h, w_h))
const I0 = CartesianIndex((r + 1, r + 1))

"""
    z_cylinder(α)

A cylinder with its axis in the xy-plane,
axis tilted α to the x-axis.
Principal curvature is (0, -r) everywhere
on the cylinder. (0,0) outside the cylinder, (∞, 0) 
on edges.  
"""
function z_cylinder(α)
    map(R) do I
        y, x = (I - I0).I
        x1 = cos(α) * x - sin(α) * y
        y1 = cos(α) * y + sin(α) * x
        sqrt(max(0, r^2 - y1^2))
    end;
end
"""
    z_cylinder_offset(α)

A cylinder with its axis in the xy-plane,
axis tilted α to the x-axis, centre at a corner.
Principal curvature is (0, -r) everywhere
on the cylinder. (0,0) outside the cylinder, (∞, 0) 
on edges.  
"""
function z_cylinder_offset(α)
    map(R) do I
        y, x = I.I
        x1 = cos(α) * x - sin(α) * y
        y1 = cos(α) * y + sin(α) * x
        sqrt(max(0, r^2 - y1^2))
    end;
end

"""
    z_sphere()

A sphere with its centre in the xy-plane.
Principal curvature is (-r, -r) everywhere
on the sphere. (0,0) outside the circle, (∞, 0) 
on edges.  
"""
function z_sphere()
    map(R) do I
        y, x = (I - I0).I
        sqrt(max(0, r^2 - x^2 - y^2))
    end
end

"""
    z_ellipsoid(; a = 0.45, tilt = 0.0)

Height field `z(x,y)`` of an ellipsoid centred in the 'xy'‑plane, with its
major half‑axis along 'y'. 

The single parameter `a ∈ [0, 0.5>` flattens the body additively:

    ay = r                # unchanged
    ax = r − a·r          # reduced once
    az = r − 2a·r         # reduced twice

`a = 0` gives a perfect sphere of radius `r`.

Returns the surface height `z ≥ 0` inside the projected ellipse and `0`
outside.

The analytical expressions for curvatures are hard to implement, so this is not a good 
choice for optimization of curvature estimation.
"""
function z_ellipsoid(; a = 0.45, tilt = 0.0)
    @assert 0 ≤ a < 0.5  "a must lie in [0, 0.5> so that all half‑axes remain positive"
    ay = r
    ax = r - a * r          # x‑axis half‑length
    az = r - 2a * r         # z‑axis half‑length
    map(R) do I
        y0, x0 = (I - I0).I
        x = y0 * sin(tilt) - x0 * cos(tilt)
        y = y0 * cos(tilt) + x0 * sin(tilt)
        s = (y^2)/(ay^2) + (x^2)/(ax^2)
        s ≥ 1 ? 0.0 : az * sqrt(1 - s)
    end
end


"""
    z_paraboloid(; a = -TestMatrices.r, b = 0.5 * TestMatrices.r)
    
    z(x, y) = (x^2)/(2a) + (y^2)/(2b)

Elliptic paraboloid, convex everywhere:  a > 0, b > 0
Hyperbolic paraboloid: Default. a and b of opposite sign. Saddle at the orgin.

Principal curvatures at (x, y) = (0, 0):

- κ₁ = 1 / a      (along x-axis)
- κ₂ = 1 / b      (along y-axis)
"""
function z_paraboloid(; a = -TestMatrices.r, b = 0.5 * TestMatrices.r)
    map(R) do I
        y, x = (I - I0).I
        (x^2)/(2a) + (y^2)/(2b)
    end
end

"""
    z_cos(; λ = r, , mult = r)    
"""
function z_cos(; λ = r, mult = r)
    map(R) do I
        y, x = (I - I0).I
        mult * cos(π * x / λ) + mult * cos(π * y / λ)
    end
end


"""
    z_exp3(; λ = r, , mult = r^2)    
"""
function z_exp3(; λ = r, mult = r)
    map(R) do I
        y, x = (I - I0).I
        mult * (x / λ)^3 + 
            mult * (y / λ)^3 
    end
end

"""
    z_plane(; λ = r, , mult = r^2)    
"""
function z_plane(; a = 0.1, b = 0.2)
    map(R) do I
        y, x = (I - I0).I
        a * x + b * y
    end
end


"""
    z_ridge_peak_valleys(;
        mult = r / 8,
        λ_arc = 0.3182291666666667,
        λ_rad = 1.875,
        spiral = 0.5)
    )
"""
function z_ridge_peak_valleys(;
    mult = r / 8,
    λ_arc = 0.3182291666666667,   # Wave length fraction of 2π inner of spiral
    λ_rad = 1.875,    # Wave length along radial in terms of half-screen width
    spiral = 0.5
    )
    map(R) do I
        negy, x = (I - I0).I # Centre origin
        y = - negy
        ρ = sqrt(float(x^2 + y^2))
        θ = mod2pi(atan(y, x))
        # Non-linear radial, normalized by half-screen-diagonal
        rnl = max(0.0,  ( ρ / √2r) - λ_rad / 6)
        # Tangential phas
        ωt = 2π * (θ * λ_arc  + spiral * rnl)
        # Phase, radial 
        ωr = ρ * 2π / (λ_rad * r)
        # Amount of cos or zigzag 
        shapepicker = (y + r)/(2r) # 1 at top to 0 at bottom
        # Interpolate between cos and zigzag
        z_norm = shapepicker * zigzag(ωt) * zigzag(ωr) + (1 - shapepicker) * cos(ωt) * cos(ωr)
        mult * z_norm
    end
end


"""
    z_wavy(; a = 0.01r, b = 0.15r, c = 0.1)
"""
function z_wavy(; a = 0.01r, b = 0.15r, c = 0.1)
    @assert ! (a == 0 &&  b == 0)
    map(R) do I
        y, x = (I - I0).I
        θ = atan(y, x)
        ρ = hypot(x, y)
        # Tangential
        zθ = sin(2θ)
        # Radial
        rwav = c * sin(8θ)
        phase = π * 2ρ / (r * (1 + rwav))
        repeatphase = mod(phase, π / 2)
        brokenphase = repeatphase < π / 4 ? repeatphase :  repeatphase + 6π / 4
        zr = -cos(brokenphase)
        θrscaledown = ρ < r / 3 ? 3ρ / r : 1.0 
        (a * zθ + b * zr) * θrscaledown 
    end
end

"""
    z_sphere_with_bulge(; a = 2096.0, b = a / 250, c = 12.0)

Centre at upper left. The default bulge (amplitude controlled by `b`) is
hardly visible. The entire output matrix is concave, but one component
varies, leading to switching of major-minor directions.
"""
function z_sphere_with_bulge(; a = 2096.0, b = a / 250, c = 12.0)
    functor = Fzx(; a, b, c)
    map(R) do I
        y, x = I.I
        rxy = hypot(x, y)
        functor(rxy)
    end
end

###############################
# Support for `z_...` functions
###############################

# For `z_ridge_peak_valleys`
function zigzag(x)
    sel = mod(x, 2π) 
    y = 2 * mod(x, π) / π - 1
    if sel ≈ 0 # Floating points can be tricky, hence special case
        1.0
    elseif sel ≈ π
        -1.0
    elseif sel < π
        - y
    else
        y
    end
end

# For `z_sphere_with_bulge`


# Find radius ρ given θ
fρθ(θ, a, b, c) = (a + b * cos(c * θ))

# Find x (i.e. rxy ) given θ
fxθ(θ, a, b, c) = fρθ(θ, a, b, c) * cos(θ)

# Find elevation z given θ"
fzθ(θ, a, b, c) = fρθ(θ, a, b, c) * sin(θ)

struct Fzx{P}
    a::Float64
    b::Float64
    c::Float64
    prob::P
end
# Constructor
function Fzx(; a = 2096.0, b = a / 250, c = 12.0, 
    θ_lower = π / 4, θ_upper = π / 2 - eps(Float64))
    xmin = fxθ(θ_upper, a, b, c)
    xmax = fxθ(θ_lower, a, b, c)
    # Residual. Parameter 'p' will be x when solving.
    f(θ, p) = clamp(first(p), xmin, xmax) - fxθ(θ, a, b, c)
    # We're adding an input 'p = x' value, 200.0, for initialization
    prob = IntervalNonlinearProblem(f, (θ_lower, θ_upper), [200.0])
    Fzx(a, b, c, prob)
end
# Callable
function (fzx::Fzx)(x)
    fzx.prob.p[1] = x
    sol = solve(fzx.prob)
    if ! successful_retcode(sol)
        throw(ErrorException("$(sol.retcode)"))
    end
    θ = Float64(sol.u)
    fzθ(θ, fzx.a, fzx.b, fzx.c)
end

##################
# Related utilties
##################

"""
    principal_curvatures_paraboloid(x, y; a = 1.0, b = 0.5a)
    principal_curvatures_paraboloid(pt::CartesianIndex; a = 1.0, b = 0.5a)

Returns the principal curvatures (κ₁, κ₂) of the elliptic paraboloid

    z(x,y) = (x²)/(2a) + (y²)/(2b)

where -1 <= x <= 1 
      -1 <= y <= 1

evaluated at point (x,y).
"""
function principal_curvatures_paraboloid(x, y; a = 1.0, b = 0.5a)
    fx = x / a
    fy = y / b
    W2 = 1 + fx^2 + fy^2
    W32 = W2^(3/2)
    H = ((1 + fy^2) * (1/a) + (1 + fx^2) * (1/b)) / (2 * W32)
    K = (1 / (a * b)) / W2^2
    sqrt_discriminant = sqrt(H^2 - K)
    κ₁ = H + sqrt_discriminant
    κ₂ = H - sqrt_discriminant
    sort([κ₁, κ₂], rev = true)
end
function principal_curvatures_paraboloid(pt::CartesianIndex; a = 1.0, b = 0.5a)
    y, x = (pt - I0).I
    principal_curvatures_paraboloid(x, y; a, b)
end



"""
    background(z::AbstractMatrix{<:Real}; Δc = -(-(extrema(z)...)) / 10, wc = Δc / 10)
    background(f::T; kws...) where T<: Union{AbstractXYFunctor,
            AbstractIJFunctor, BidirectionAtXY, BidirectionInDomain,
            Vec2InDomain}

An image of the same size as input. Argument `z` is shown as an elevation matrix with contour bands.

# Keyword arguments

- `Δc` Contour lines distance
- `wc` Width of contour lines along the z axis.
"""
function background(z::AbstractMatrix{<:Real}; Δc = -(-(extrema(z)...)) / 10, wc = Δc / 10)
    foo = scaleminmax(extrema(z)...)
    scaled_z = foo.(z)
    img = map(scaled_z) do z
        RGB{N0f8}(get(PALETTE_BACKGROUND, z))
    end
    # Add simple contour lines, too
    map!(img, z, img) do zz, pix
        mod(zz, Δc) < wc ? RGB{N0f8}(0.1, 0.1, 0.1) : pix 
    end
end
function background(f::T; kws...) where T<: Union{AbstractXYFunctor,
            AbstractIJFunctor, BidirectionAtXY, BidirectionInDomain,
            Vec2InDomain}
    background(z_matrix(f); kws...)
end
end # module


