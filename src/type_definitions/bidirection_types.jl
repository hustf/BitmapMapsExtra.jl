 # This defines functor types.
 # BidirectionOnGrid
 # BidirectionInDomain
 # BidirectionAtXY
 # SelectedVec2AtXY



struct BidirectionOnGrid{F, T, LC} <: AbstractIJFunctor
    fdir!::F # e.g 𝐊!. Input signature: (K, vα, vβ, vκ, P, M, vϕ)
    z::Matrix{T} # Image-sized
    Ω::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}
    vα::MVector{4, Float64}
    vβ::MVector{2, Float64}
    vκ::MVector{4, Float64}
    P::MMatrix{3, 3, Float64, 9}
    lpc::LC
    K::MMatrix{2, 2, Float64, 4}
end
# Constructor
function BidirectionOnGrid(fdir!, z)
    _, Ω, _, P, K, vα, vκ, vβ, lpc = allocations_curvature(CartesianIndices(z))
    BidirectionOnGrid(fdir!, z, Ω, vα, vβ, vκ, P, lpc, K)
end
# Callable, returns a tensormap K for the pixel or cell
function (b::BidirectionOnGrid)(i::Int64, j::Int64)
    rngi = b.Ω.indices[1] .+ i
    rngj = b.Ω.indices[2] .+ j
    vi = view(b.z, rngi, rngj)
    b.fdir!(b.K, b.vα, b.vβ, b.vκ, b.P, vi, VΦ, b.lpc)
end
@define_show_with_fieldnames BidirectionOnGrid

# Intermediate level
struct BidirectionInDomain{F, T, LC}
    bdog::BidirectionOnGrid{F, T, LC}
    corners::SMatrix{2, 2, TENSORMAP, 4}
    normalize::Bool
    minnorm::Float64
    lastvalue::TENSORMAP
end
# Constructor
function BidirectionInDomain(fdir!, z::Matrix{<:AbstractFloat}; normalize = false, minnorm = 1e-4)
    @assert minnorm >= 0
    # Can generate a zero corner value (a 4x4), mutable
    # Establish four mutable corners in a static matrix 
    corners = SMatrix{2, 2, TENSORMAP, 4}([TENSORMAP(zeros(Float64, 2, 2)) for i in 1:2, j in 1:2] )
    lastvalue = TENSORMAP(zeros(Float64, 2, 2))
    # Store as BidirectionInDomain
    BidirectionInDomain(BidirectionOnGrid(fdir!, z), corners, normalize, minnorm, lastvalue)
end
# Callable, returns K, a TENSORMAP for floating point coordinates
function (bid::BidirectionInDomain)(jfloat, ifloat)
    # Identify (up to four) points on grid
    j1, j2 =  Int(floor(jfloat)), Int(ceil(jfloat))
    i1, i2 =  Int(floor(ifloat)), Int(ceil(ifloat))
    di = ifloat - i1
    dj = jfloat - j1
    update_corners!(bid, i1, j1, i2, j2, di, dj, bid.minnorm)
    # Note, the dj, di are in "(x, y)" order.
    bid.lastvalue .= interpolate_and_normalize_directions!(bid.lastvalue, bid.corners, dj, di, bid.normalize, bid.minnorm)
    bid.lastvalue
end
@define_show_with_fieldnames BidirectionInDomain

"""
    struct BidirectionAtXY{F, T, LC}
        bid::BidirectionInDomain{F, T, LC}
        negy::NegateY
        d::Domain
        K::TENSORMAP
        lastvalue::MVector{2, Float64}
        major::Ref{Bool}
        major_start::Bool
    end

Calling an instance with coordinates will return a 180° symmetric vector.
"""
struct BidirectionAtXY{F, T, LC}
    bid::BidirectionInDomain{F, T, LC}
    negy::NegateY
    d::Domain
    K::TENSORMAP
    lastvalue::MVector{2, Float64}
    major::Ref{Bool} # Perhaps not the best solution.
    major_start::Bool
end
# Constructor
function BidirectionAtXY(fdir!, z, major::Bool; kws...)
    bid = BidirectionInDomain(fdir!, z; kws...)
    negy = NegateY(z)
    R = CartesianIndices(z)
    Ω = CartesianIndices((-2:2, -2:2))
    d = Domain(R, Ω)
    K = TENSORMAP(zeros(Float64, 2, 2))
    lastvalue = K[:, 1]
    BidirectionAtXY(bid, negy, d, K, lastvalue, Ref{Bool}(major), major)
end
# Callable. The returned 2d vector is ambiguous, i.e. "180°-symmetric". 
function (baxy::BidirectionAtXY)(x::Float64, y::Float64)
    if baxy.d(x, y)
        baxy.K .= baxy.bid(x, baxy.negy(y))
        if baxy.major[]
            baxy.lastvalue .= baxy.K[:, 1]
        else
            baxy.lastvalue .= baxy.K[:, 2]
        end
    else
        # (x, y) is out of domain. Solvers make such calls fairly often.
        # We don't want an early fail because of that.
        baxy.K .= 0.0
        baxy.lastvalue .= baxy.K[:, 1]
    end
    baxy.lastvalue
end
@define_show_with_fieldnames BidirectionAtXY




"""
    struct SelectedVec2AtXY{F, T, LC} <: AbstractXYFunctor
        baxy::BidirectionAtXY{F, T, LC}
        v::MVector{2, Float64}
        negy::NegateY # Refers baxy.negy
        flip::Ref{Bool}
        flip_start::Bool
    end

Functor for making streamlines from "tensor maps", namely
two 180° symmetrical 2d vectors. One such bidirectional 
2d vector is 'major', one is `minor`.

From any single point, we could integrate along either 
of four directions.

# Constructor

Use the constructor to pick one direction:

    saxy = SelectedVec2AtXY(fdir!, z, major::Bool, flip::Bool)

# Arguments

- `fdir!` is your function, which must have the same function signature as `𝐊!`. 
  Most of the arguments are just buffers which you may use or not.
- `z` is a matrix of floats, such as an elevation map.
- `major` = true selects the direction with the signed largest value.
- `flip` = false will integrate along this curve to the positive side.

While integrating along one such curve, the value may change sign. We would want to
continue without re-tracing previously visited point. Hence, if a sign change is detected,
the `flip` value itself will also flip. The value when constructed can be restored by 
calling `reset!(saxy)`. When plotting multiple streamlines, this functor is reset at
every starting point.

# Callable

    saxy(x, y)
    --> MVector{2, Float64}

...where x and y are floating point coordinates in an x-y up coordinate system.
"""
struct SelectedVec2AtXY{F, T, LC} <: AbstractXYFunctor
    baxy::BidirectionAtXY{F, T, LC}
    v::MVector{2, Float64}
    negy::NegateY # Refers baxy.negy
    flip::Ref{Bool}
    flip_start::Bool
end
# Constructor
function SelectedVec2AtXY(fdir!, z, major::Bool, flip::Bool; normalize = true, minnorm = 1e-4)
    v = MVector{2, Float64}([0.0, 0.0])
    baxy = BidirectionAtXY(fdir!, z, major; normalize, minnorm)
    negy = baxy.negy
    SelectedVec2AtXY(baxy, v, negy, Ref{Bool}(flip), flip)
end
# Callable. The returned 2d vector is "360°-symmetric"
function (saxy::SelectedVec2AtXY)(x::Float64, y::Float64)
    if saxy.flip[]
        saxy.v .=- saxy.baxy(x, y)
    else
        saxy.v .= saxy.baxy(x, y)
    end
    saxy.v
end
@define_show_with_fieldnames SelectedVec2AtXY


# Method extensions
z_matrix(fij::BidirectionOnGrid) = fij.z
z_matrix(fxy::BidirectionInDomain) = z_matrix(fxy.bdog)
z_matrix(fxy::BidirectionAtXY) = z_matrix(fxy.bid)
z_matrix(saxy::SelectedVec2AtXY) = saxy.baxy.bid.bdog.z

function size(f::T) where T<: Union{BidirectionAtXY, BidirectionInDomain}
    size(z_matrix(f))
end



# Resetting will not destroy the inherent z-data,
# just prepare for integrating from a new starting point.
function reset!(saxy::SelectedVec2AtXY)
    saxy.flip[] = saxy.flip_start    
    saxy.baxy.major[] = saxy.baxy.major_start
    saxy.v .= 0.0
    saxy.baxy.K .= 0.0
    saxy
end
reset!(fxy::AbstractXYFunctor) = fxy