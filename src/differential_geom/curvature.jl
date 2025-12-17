# Contains `principal_curvature_components` and it's non-allocating cousin, 
# as well as supporting functions.
# Also contains `allocations_curvature`.
# Relies on `tangent_basis.jl` and `bidirectional_quantity.jl`.
# Relies on constants KERN´´ and VΦ.

"""
    principal_curvature_components(z::Matrix{<:AbstractFloat}, pt::CartesianIndex)

Note that there is an alias 𝐊 for brevity.
Use the '!' methods where speed matters.
"""
function principal_curvature_components(z::Matrix{<:AbstractFloat}, pt::CartesianIndex)
    Ri, Ω, v, P, K, vα, vκ, vβ, lpc = allocations_curvature(CartesianIndices(z))
    if ! (pt ∈ Ri)
        # Too close to edge of z    
        K .= NaN
    else
        # Window 
        win = view(z, Ω .+ pt)
        # Update K etc. 
        principal_curvature_components!(K, vα, vβ, vκ, P, win, VΦ, lpc)
    end
    K
end

"""
    principal_curvature_components!(K, vα, vβ, vκ, P, M, vϕ, lpc)

Alias: 𝐊!.

- `K` is the result value (a 2x2 TENSORMAP), updated in place.
- `M` is a 5x5 window centered on the current evaluation element. 
- `vϕ` are sample angles, normally the optimized values VΦ.
- `vα`, `vβ`, `vκ` and `P` are buffers for intermediate calculations. 
  These will also be mutated. 
  
Use `allocations_curvature` for generating all the arguments in one go.
"""
function principal_curvature_components!(K, vα, vβ, vκ, P, M, vϕ, lpc)
    _, κ1, κ2 = principal_curvatures_and_angles!(vβ, vα, vκ, P, M, vϕ, lpc)
    # Express curvatures κ and directions β as 4x4 matrix K.
    # K is a second-order tensor's components in a screen-aligned basis.
    components_matrix!(K, κ1, κ2, vβ)
end

"""
    principal_curvature_normalized!(Kᵤ, vα, vβ, vκ, P, M, vϕ, lpc; minnorm = 1e-4)

Alias: Kᵤ!.

See `principal_curvature_components!` 

Use this for making `BidirectionOnGrid` objects for plotting tensormap glyphs.

NOTE: Avoid this for making `BidirectionInDomain`, `BidirectionAtXY` and `SelectedVec2AtXY` 
objects. Instead, use keywords `normalize = true`, possibly combined with e.g. `minnorm = 1e-4`

Both directions are normalized to length 1, or set to zero when below threshold `minnorm`. 
The sign remains unchanged.
"""
function principal_curvature_normalized!(Kᵤ, vα, vβ, vκ, P, M, vϕ, lpc; minnorm = 1e-4)
    principal_curvature_components!(Kᵤ, vα, vβ, vκ, P, M, vϕ, lpc)
    normalize_or_zero!(Kᵤ, minnorm)
end

"""
    principal_curvatures_and_angles!(vβ, vα, vκ, P, M, vϕ, lpc)
    --> vβ, κ1, κ2

Internal function.
vβ is modified in-place. So is vα, vκ, P and lpc. A bit of a mess, really.
"""
function principal_curvatures_and_angles!(vβ, vα, vκ, P, M, vϕ, lpc)
    @assert size(vα) == (4,)
    @assert size(vβ) == (2,)
    @assert size(vκ) == (4,)
    @assert size(P) == (3, 3)
    @assert size(M) == (5, 5)
    # Find P in-place. vβ is a temporary storage here. 
    tangent_basis!(P, vβ, M)
    # vϕ gives the sampling direction in the tangent plane
    # described by P. Sampling must be done in the screen plane,
    # since our data is an elevation matrix with indices as -y, x
    # 
    # Find 'screen angles' vα which are projections of the fixed vϕ
    angle_tangent_to_xy!(vα, vϕ, P)
    # Sample curvatures vκ
    sample_curvature_in_directions!(vκ, M, P, vα)
    # Calculations in the tangent plane
    κ1, κ2, ϕ1 = principal_curvature_and_direction(vκ, lpc)
    @assert typeof(κ1) <: Float64
    @assert typeof(κ2) <: Float64
    @assert typeof(ϕ1) <: Float64
    # Principal angles vβ in the yx (screen) plane. 
    tangent_dyad_to_xy!(vβ, ϕ1, P)
    vβ, κ1, κ2
end

"""
    sample_curvature_in_directions!(vκ, M, P, vα)
"""
function sample_curvature_in_directions!(vκ, M, P, vα)
    # Tangent plane normal's components:
    𝐧ₜₚ₁ = P[7]
    𝐧ₜₚ₂ = P[8]
    𝐧ₜₚ₃ = P[9]
    # This gives us an estimate for the derivative of the surface 
    # in the plane intersected by α
    # ζ´ = -(cos(α) * 𝐧ₜₚ₁ + sin(α) * 𝐧ₜₚ₂) / 𝐧ₜₚ₃
    # Mutates vκ according to vα
    vκ .= map(α -> sample_in_direction(α, M, 𝐧ₜₚ₁, 𝐧ₜₚ₂, 𝐧ₜₚ₃), vα)
    vκ
end

function sample_in_direction(α, M, 𝐧ₜₚ₁, 𝐧ₜₚ₂, 𝐧ₜₚ₃)
    # Shorthand
    co = cos(α)
    si = sin(α)
    #
    # Normal curvature and unit tangent vector in the 'matrix coordinate system',
    # where when α = 0 : 
    #           local axis lx points along increasing rows in the unrotated matrix (that is: down in an image).
    #           Sampling from increasing rows could be faster than sampling from decreasing rows because of how
    #           memory access works.
    #       when α = π / 4: 
    #           local axis lx points along decreasing columns (that is: left in an image)
    # First derivative, taken from the tangent plane
    z´ = -(co * 𝐧ₜₚ₁ + si * 𝐧ₜₚ₂) / 𝐧ₜₚ₃
    # Second derivative
    z´´ = sampled_second_derivative_in_direction(M, co, -si)
    # Extrinsic curvature in this direction (from Kreyszig AEM)
    z´´ / (1 + z´^2)^(3 / 2)
end

function sampled_second_derivative_in_direction(M, co, si)
    @assert size(M) == (5, 5)
    accum = 0.0
    for pos in 1:5
        k = KERN´´[pos]
        x = 3 + (pos - 3) * co
        y = 3 + (pos - 3) * si
        accum += k * sample_at_float(M, x, y)
    end
    accum
end

function sample_at_float(M, x, y)
    # Identify (up to four) points on grid around (x, y)
    j1, j2 =  Int(floor(x)), Int(ceil(x))
    i1, i2 =  Int(floor(y)), Int(ceil(y))
    # Coordinates within unit square
    u = x - j1
    v = y - i1
    # Interpolate
    (1 - u) * (1 - v) * M[i1, j1] + 
          (1 - u) * v * M[i2, j1] +
          u * (1 - v) * M[i1, j2] + 
                u * v * M[i2, j2]
end

"""
    angle_xy_to_tangent!(vϕ, vα, P)

This works with tangent bases with the condition that
the tangent plane's primary axis lies in the yz plane, and 
the tangent plane's third axis is flipped in the direction of z.

Projects angles from the x-axis around the z-axis into the tangent
plane, along the z-axis. The projected angle similarly refers to the 
tangent plane's primary axis.

Reversible transformation, see `angle_tangent_to_xy!`
"""
function angle_xy_to_tangent!(vϕ, vα, P)
    a, _, c, d, e, f, g, h, i = P
    @assert i !== 0
    map!(vϕ, vα) do α
        co = cos(α)
        si = sin(α)
        z = -(g·co + h·si) / i
        u = a·co +        c⋅z   # component along 𝐞₁
        v = d·co + e·si + f·z  # component along 𝐞₂
        ϕ = atan(v, u)
    end
end

"""
    angle_tangent_to_xy!(vα, vϕ, P)

This works with tangent bases with the condition that
the tangent plane's primary axis lies in the yz plane, and 
the tangent plane's third axis is flipped in the direction of z.

Projects angles from the tangent plane's primary axis around its third 
axis into the xy plane, along the z-axis. The projected angle similarly
refers to the x axis.

Reversible transformation, see `angle_xy_to_tangent!`
"""
function angle_tangent_to_xy!(vα, vϕ, P)
    # This works with tangent bases with 
    # the conditions we set on orientation.
    map!(vα, vϕ) do ϕ
        # Tangent-space vector [u, v, w]
        u = cos(ϕ)
        v = sin(ϕ)
        # Standard-basis vector 𝐩 = P ⋅ 𝐮
        x = P[1]·u + P[4]·v
        y = P[5]·v
        α = atan(y, x)
    end
end

"""
    tangent_dyad_to_xy!(vα, ϕ::AbstractFloat, P)

Saves a memory allocation compared to the equivalent

```
angle_tangent_to_xy!(vα, [ϕ, ϕ + π / 2], P)
```
"""
function tangent_dyad_to_xy!(vα, ϕ::AbstractFloat, P)
    # Tangent plane direction
    u = cos(ϕ)
    v = sin(ϕ)
    # Standard-basis vector 𝐩 = P ⋅ 𝐮
    # Project onto standard basis xy.
    # This works with tangent bases with 
    # the conditions we set on orientation.
    x = P[1]·u + P[4]·v
    y = P[5]·v
    vα[1] = atan(y, x)
    # In the tangent plane, the minor principal direction
    # is perpendicular to ϕ. The unit directional vectors 
    # form an orthonormal dyad.
    u, v = -v, u
    x = P[1]·u + P[4]·v
    y = P[5]·v
    vα[2] = atan(y, x)
    vα
end


"""
    principal_curvature_and_direction(vκ::T, lpc) where T <: SVector{4, Float64}
    principal_curvature_and_direction(vκ, lpc) 
    -> (Float64, Float64, Float64)

Compute the principal curvatures κ₁ ≥ κ₂ and the principal direction angle φₚ from three curvature samples.

Given measurements κ(φᵢ) at angles φᵢ (in a common tangent-plane frame), we model

    κ(φ) = κ₁ cos(φ - φₚ)^2 + κ₂ sin(φ - φₚ)^2

where φₚ is the angle (in the same frame) at which the maximum curvature κ₁ occurs.

Expanding via trigonometric identities:
    cos^2(φ - φₚ) = [1 + cos(2φ - 2φₚ)]/2
    sin^2(φ - φₚ) = [1 - cos(2φ - 2φₚ)]/2

leads to the equivalent form

    κ(φ) = a + b cos(2φ) + c sin(2φ)
    where
      a = (κ₁ + κ₂) / 2
      r = (κ₁ - κ₂) / 2  (so b = r cos(2φₚ), c = r sin(2φₚ))

We solve the linear system

    A * [a b c]ᵀ = vκ,

where

    A = [1  cos(2φ₁)  sin(2φ₁);
         1  cos(2φ₂)  sin(2φ₂);
         1  cos(2φ₃)  sin(2φ₃)
         1  cos(2φ₄)  sin(2φ₄)]
      = hcat(ones(Float64, 4), cos.(2vφ), sin.(2vφ))

which we solve with an optimized equivalent to:

    a, b, c = A \\ vκ

Then:
    r = √(b² + c²)
    κ₁ = a + r
    κ₂ = a - r
    φₚ = ½ ⋅ atan2(c, b)

Ensuring φₚ is in the correct branch of [0, 2π).

# Parameters
- vκ    Three curvature samples κ(φᵢ)
- vφ    Three angles φᵢ in the same tangent-plane frame
- lpc   Linear problem reuseable cache

# Returns
- κ₁    Maximum principal curvature
- κ₂    Minimum principal curvature
- φₚ    Principal-direction angle (same units/frame as vφ)


# References
1. https://en.wikipedia.org/wiki/Euler%27s_theorem_(differential_geometry)
"""
function principal_curvature_and_direction(vκ::T, lpc) where T <: SVector{4, Float64}
    # Update and solve A * [a; b; c] = vκ
    lpc.b = vκ
    a, b, c = solve!(lpc)
    # Principal‐curvature values from a, b, c
    r = hypot(b, c)       # faster than √(b^2 + c^2)
    # TEMP 
    if r < 1*10^-4.2
        a = 0.0
        b = 0.0
        c = 0.0
        r = 0.0
    end
    κ1 = a + r            # maximum principal curvature
    κ2 = a - r            # minimum principal curvature
    # First maximum of I) equals κ1 and occurs at (from differentiation)
    ϕp = atan(c, b) / 2
    if isnan(ϕp)
        throw(ErrorException("Unexpected NaN: vκ, a, b, c = $vκ $a $b $c"))
    end
    @assert κ1  >= κ2
    κ1, κ2, ϕp
end
principal_curvature_and_direction(vκ, lpc) = 
    principal_curvature_and_direction(SVector{4, Float64}(vκ), lpc)











"""
    allocations_curvature(R::CartesianIndices)

Allocate once, re-use at every point. This function is re-used in calculations 
which require fewer such containers than curvature.

# Arguments

R             e.g. CartesianIndices(size(z))
directions    e.g. [] or [2] or [1, 2]

# Output

Ri, Ω, v, P, K, vα, vκ, vβ, lpc

where

typeof(Ri) = CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}
typeof(Ω) = CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}
typeof(v) = StaticArraysCore.MVector{2, Float64}
typeof(P) = StaticArraysCore.MMatrix{3, 3, Float64, 9}
typeof(K) = StaticArraysCore.MMatrix{2, 2, Float64, 4}
typeof(vα) = StaticArraysCore.MVector{4, Float64}
typeof(vκ) = StaticArraysCore.MVector{4, Float64}
typeof(vβ) = StaticArraysCore.MVector{2, Float64}
typeof(lpc) = LinearCache{...}
"""
function allocations_curvature(R::CartesianIndices; vϕ = VΦ)
    # Define an internal domain Ri, since simple padding options
    # would not yield interesting curvature anyway.
    n = 5
    side = n ÷ 2
    is, js = side + 1, side + 1
    ie, je = size(R,1) - side - 1, size(R,2) - side - 1
    Ri = CartesianIndices((is:ie, js:je))
    # Relative indices, for defining a sliding window around each point
    Ω = CartesianIndices((-side:side, -side:side))
    # Pre-allocate mutable containers for re-using memory. 
    #   fixed-sized 2d vector, multiple purpose
    v = MVector{2, Float64}(Array{Float64, 1}(undef, 2))
    #   statically-sized, tangent basis
    P = MMatrix{3, 3, Float64}(Array{Float64, 2}(undef, (3, 3)))
    #  fixed-sized, curvature components 
    K = MMatrix{2, 2, Float64}(Array{Float64, 2}(undef, (2, 2)))
    #  fixed-sized, curvature sample angles in screen plane 
    vα = MVector{4, Float64}(Array{Float64, 1}(undef, 4))
    #  fixed-sized, sample curvature values  
    vκ = MVector{4, Float64}(Array{Float64, 1}(undef, 4))
    # Principal and minor principal angles. 
    vβ = MVector{2, Float64}(Array{Float64, 1}(undef, 2))
    # For finding which principal directions and size could 
    # lead to our samples. We prepare by constructing a design matrix A
    # and bake it into a cache for solving A * [a; b; c] = vκ.
    # 
    # See `principal_curvature_and_direction` regarding the equations.
    A = SMatrix{4, 3, Float64, 12}(hcat(ones(Float64, size(vϕ, 1)), cos.(2vϕ), sin.(2vϕ)))
    vκs = SVector{4, Float64}(Array{Float64, 1}(undef, 4))
    lpc = init(LinearProblem{true}(A, vκs), NormalCholeskyFactorization(), OperatorAssumptions(false, condition=OperatorCondition.WellConditioned))
    #
    Ri, Ω, v, P, K, vα, vκ, vβ, lpc
end






