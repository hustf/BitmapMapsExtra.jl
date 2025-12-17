 # This builds on domain, is a requisite for
 # streamlines.jl,  and define types for
 # 'non-discretization' of matrices, making them
 # work like continuous functions.
 #
 # Defines types:

 # Vec2AtXY (top level)
 # Vec2OnGrid
 # Vec2InDomain

 """
Supertype for functors with an on grid, discrete i, j domain,
where i is pixel index from image top, j index is image right.
Used for dispatch in plotting.
"""
abstract type AbstractIJFunctor end


struct Vec2OnGrid{F, T} <: AbstractIJFunctor
    Ω::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}# Square neighborhood
    fdir!::F # e.g 𝐧ₚ!(v, M)
    z::Matrix{T} # Image-sized
    lastvalue::MVector{2, Float64}
end
# Constructor
function Vec2OnGrid(fdir!, z)
    Ω = CartesianIndices((-2:2, -2:2))
    lastvalue = zero(MVector{2, Float64})
    Vec2OnGrid(Ω, fdir!, z, lastvalue)
end
# Callable, returns a 2d vector
function (vog::Vec2OnGrid)(i::Int64, j::Int64)
    rngi = vog.Ω.indices[1] .+ i
    rngj = vog.Ω.indices[2] .+ j
    vi = view(vog.z, rngi, rngj)
    vog.fdir!(vog.lastvalue, vi)
    vog.lastvalue
end
@define_show_with_fieldnames Vec2OnGrid


struct Vec2InDomain{F, T}
    vog::Vec2OnGrid{F, T}
    lastvalue::MVector{2, Float64}
end
# Constructor
function Vec2InDomain(fdir!, z::Matrix{<:AbstractFloat})
    vog = Vec2OnGrid(fdir!, z)
    lastvalue = zero(MVector{2, Float64})
    Vec2InDomain(vog, lastvalue)
end
# Callable, returns a 2d vector. Function of in-domain float coordinates. 'negy' is 0 at bottom of image.
function (did::Vec2InDomain)(x, negy)
    # Identify (up to four) points on grid
    j1, j2 =  Int(floor(x)), Int(ceil(x))
    i1, i2 =  Int(floor(negy)), Int(ceil(negy))
    # Shorthand 
    vog = did.vog
    lastvalue = did.lastvalue
    #
    if i1 != i2 && j1 != j2
        # Both pairs are different
        lastvalue .=  corner_weight(1, 1, x - j1, negy - i1) .* vog(i1, j1)
        lastvalue .+= corner_weight(2, 1, x - j1, negy - i1) .* vog(i2, j1)
        lastvalue .+= corner_weight(1, 2, x - j1, negy - i1) .* vog(i1, j2)
        lastvalue .+= corner_weight(2, 2, x - j1, negy - i1) .* vog(i2, j2)
    elseif i1 == i2 && j1 != j2
        # i1 equals i2, but j1 differs from j2
        lastvalue .=  linear_weight(1, x - j1) .* vog(i1, j1)
        lastvalue .+= linear_weight(2, x - j1) .* vog(i1, j2)
    elseif j1 == j2 && i1 != i2
        lastvalue .=  linear_weight(1, negy - i1) .* vog(i1, j1)
        lastvalue .+= linear_weight(2, negy - i1) .* vog(i2, j1)
        # j1 equals j2 
    else
        # i1 equals i2 and j1 equals j2
         lastvalue .= vog(i1, j1)
    end
    lastvalue
end
@define_show_with_fieldnames Vec2InDomain

"""
    abstract type AbstractXYFunctor end

Supertype for functors with a continuous x-y up domain space.
Used for dispatch in differential equation solving.

See `subtypes(AbstractXYFunctor)`.
"""
abstract type AbstractXYFunctor end

"""
    struct Vec2AtXY{F, T} <: AbstractXYFunctor
        did::Vec2InDomain{F, T}
        negy::NegateY
        d::Domain
        v::MVector{2, Float64}
    end

Functor for making streamlines from 2d vector fields.

# Constructor:

    vaxy = Vec2AtXY(fdir!, z)

# Arguments

- `fdir!` is your function, which must have the same function signature as `descent!`. 
- `z` is a matrix of floats, such as an elevation map.

Callable: 

    vaxy(x, y) 
    --> MVector{2, Float64}

...where x and y are floating point coordinates in an x-y up coordinate system.
"""
struct Vec2AtXY{F, T} <: AbstractXYFunctor
    did::Vec2InDomain{F, T}
    negy::NegateY
    d::Domain
    v::MVector{2, Float64}
end
# Constructor
function Vec2AtXY(fdir!, z)
    did = Vec2InDomain(fdir!, z)
    negy = NegateY(z)
    R = CartesianIndices(z)
    Ω = CartesianIndices((-2:2, -2:2))
    d = Domain(R, Ω)
    v = MVector{2, Float64}([0, 0])
    Vec2AtXY(did, negy, d, v)
end
# Callable. y is positive up on screen. Returns [0.0,0.0]
# when out-of domain (the callback set doesn't immediately
# stop calculation).
function (vaxy::Vec2AtXY)(x, y)
    if vaxy.d(x, y)
        vaxy.v .= vaxy.did(x, vaxy.negy(y))
    else
        # This occurs fairly often, though isn't presented as part of the solution.
        #printstyled("::Vec2AtXY(x, y) == ($x, $y) ∉ $(vaxy.d) \n", color =:176)
        vaxy.v .= 0.0
        vaxy.v
    end
    vaxy.v
end
@define_show_with_fieldnames Vec2AtXY


# Method (extensions)

z_matrix(fij::Vec2OnGrid) = fij.z
z_matrix(fij::Vec2InDomain) = z_matrix(fij.vog)
z_matrix(fxy::Vec2AtXY) = z_matrix(fxy.did)

function size(f::T) where T<: Union{AbstractXYFunctor, AbstractIJFunctor, 
    Vec2InDomain}
    size(z_matrix(f))
end
