# File contains utilty functions for
# the 4x4 mutable matrix ('K') we use as a tensor map
# for the basis.  

##################################################
# Interpolation while respecting major-minor order
##################################################

"""
    interpolate_and_normalize_directions!(K::TENSORMAP, 
        corners::SMatrix{2, 2, TENSORMAP}, 
        xus::AbstractFloat, 
        yus::AbstractFloat, 
        normalize::Bool, 
        minnorm::Float64)
    --> TENSORMAP

Interpolates between four TENSORMAPs (each contaiing two bidirectional 2d vectors) at the corners of the unit square, 
and then normalizes each column to unit length or zero.

Why normalize? Already in the simpler 1D case of linear interpolation between two 2d vectors 
for t ∈ (0, 1), we have, by the strict triangle inequality,
    
    ||(1 -t)𝐯₀ + t𝐯₁|| < (1 -t)||𝐯₀|| + t||𝐯₁|| for non-parallel 𝐯₀, 𝐯₁.

Thus, interpolated vectors in the interior tend to be shorter than the corner vectors. This could easily be misread
as a loss of "strength" where corner vectors diverge. We therefore keep bilinear interpolation of direction 
but discard interpolated magnitude by normalizing to unit length.

Finally, normalization of close to zero-magnitude vectors would be equally misleading. Hence, columns with
`magnitude < minnorm` are set to zero.

# CONSIDER TODO: Separate functions for normalize true or false.
"""
function interpolate_and_normalize_directions!(K::TENSORMAP, corners::SMatrix{2, 2, TENSORMAP}, xus::AbstractFloat, yus::AbstractFloat, normalize::Bool, minnorm::Float64)
    x = clamp(xus, 0.0, 1.0)
    y = clamp(yus, 0.0, 1.0)
    if xus !== x
        throw(" ooo: x")
    end
    if yus !== y
        throw(" ooo: y")
    end    
    f11 = corners[1, 1] # (x, y) = (0, 0) 
    f12 = corners[1, 2] # (x, y) = (0, 1)
    f21 = corners[2, 1] # (x, y) = (1, 0)
    f22 = corners[2, 2] # (x, y) = (1, 1)
    K .= (1 - x) * (1 - y) .* f11 .+ 
                   (1 - x) * y .* f21 .+
                   x * (1 - y) .* f12 .+ 
                         x * y .* f22
    if normalize
        normalize_or_zero!(K, minnorm)
    else
        # Even if we don't normalize, we compensate the
        # lengths of each bidirectional vector
        @inbounds for j in axes(K, 2)
            l = norm(K[:, j])
            if l > MAG_EPS
                l11 = norm(f11[:, j])
                l12 = norm(f12[:, j])
                l21 = norm(f21[:, j])
                l22 = norm(f22[:, j])
                li = (1 - x) * (1 - y) .* l11 .+ 
                           (1 - x) * y .* l21 .+
                           x * (1 - y) .* l12 .+ 
                                x * y .* l22
                mult = li / l
                K[:, j] .*= mult
            end
        end
        enforce_minnorm!(K, minnorm)
    end 
    K
end

"""
    corner_weight(i, j, x::Float64, y::Float64)
    --> Float64

For interpolation in a unit square
"""
function corner_weight(i, j, x::Float64, y::Float64)
    i == 1 && j == 1 && return (1 - x) * (1 - y)
    i == 2 && j == 1 && return (1 - x) * y
    i == 1 && j == 2 && return x * (1 - y)
    i == 2 && j == 2 && return x * y
    throw(ErrorException("$i $j"))
end

"""
    linear_weight(i, j, x::Float64, y::Float64)
    --> Float64

For interpolation along a line 0..1.
"""
function linear_weight(j, x::Float64)
    j == 1 && return (1 - x)
    j == 2 && return x
    throw(ErrorException("$j $x"))
end
"""
    update_corners!(bid::BidirectionInDomain, i1, j1, i2, j2, di, dj)

Calculate and store tensormap K for closest on-grid coordinates.
Purpose is linear interpolation.

Major and minor principal  directions are not allowed to swap 
column order between corners. This is to ensure
that such a swap is detectable, even with interpolation.


"""
function update_corners!(bid::BidirectionInDomain, i1, j1, i2, j2, di, dj, minnorm)
    # Alias
    mc = bid.corners
    bdog = bid.bdog
    # Fallback corner when another must be dropped
    dominant_i = 1 + Int(round(di))
    dominant_j = 1 + Int(round(dj))
    if i1 != i2 && j1 != j2     # Both pairs are different
        update_corners_all!(mc, bdog, i1, j1, i2, j2)
        enforce_compatibility_with_dominant!(mc, dominant_i, dominant_j, minnorm)
    elseif i1 == i2 && j1 != j2 # i1 equals i2, but j1 differs from j2
        update_corners_i_equal!(mc, bdog, i1, j1, j2)
        enforce_compatibility_with_dominant!(mc, dominant_i, dominant_j, minnorm)
    elseif j1 == j2 && i1 != i2 # j1 equals j2
        update_corners_j_equal!(mc, bdog, i1, j1, i2)
        enforce_compatibility_with_dominant!(mc, dominant_i, dominant_j, minnorm)
    else                        # i1 equals i2 and j1 equals j2
        update_corners_ij_equal!(mc, bdog, i1, j1) 
        # No interpolation, no compatibilty problem
    end
end

function update_corners_all!(mc, bdog, i1, j1, i2, j2)
    @debug "update corners all"
    mc[1, 1] .= bdog(i1, j1)
    mc[2, 1] .= bdog(i2, j1)
    mc[1, 2] .= bdog(i1, j2)
    mc[2, 2] .= bdog(i2, j2)
    mc
end
function update_corners_i_equal!(mc, bdog, i, j1, j2)
    @debug "update corners i equal"
    mc[1, 1] .= bdog(i, j1)
    mc[2, 1] .= mc[1, 1]
    mc[1, 2] .= bdog(i, j2)
    mc[2, 2] .= mc[1, 2]
    mc
end
function update_corners_j_equal!(mc, bdog, i1, j, i2)
    @debug "update corners j equal"
    mc[1, 1] .= bdog(i1, j)
    mc[2, 1] .= bdog(i2, j)
    mc[1, 2] .= mc[1, 1]
    mc[2, 2] .= mc[2, 1]
    mc
end
function update_corners_ij_equal!(mc, bdog, i, j)
    @debug "update corners ij equal"
    mc[1, 1] .= bdog(i, j)
    mc[2, 1] .= mc[1, 1]
    mc[1, 2] .= mc[1, 1]
    mc[2, 2] .= mc[1, 1]
    mc
end


function enforce_compatibility_with_dominant!(mc, dominant_i::Integer, dominant_j::Integer, minnorm)
    enforce_zero_compatibility_with_dominant!(mc, dominant_i, dominant_j, 1, minnorm)
    enforce_zero_compatibility_with_dominant!(mc, dominant_i, dominant_j, 2, minnorm)
    enforce_column_compatibility_with_dominant!(mc, dominant_i, dominant_j)
    enforce_direction_compatibility_with_dominant!(mc, dominant_i, dominant_j, 1)
    enforce_direction_compatibility_with_dominant!(mc, dominant_i, dominant_j, 2)
    mc
end

"""
    enforce_zero_compatibility_with_dominant!(mc, dominant_i, dominant_j, column, minnorm)

If the dominant corner interpolation contains zero-vectors, set all to zero. 
If dominant is non-zero, set zero valued corners equal to the dominant value.
"""
function enforce_zero_compatibility_with_dominant!(mc, dominant_i, dominant_j, column, minnorm)
    dom = mc[dominant_i, dominant_j][:, column]
    if norm(dom) < minnorm
        # Set all non-dominant corners to zero too
        @inbounds for j in axes(mc, 2), i in axes(mc, 1)
            i == dominant_i && j == dominant_j && continue
            @debug "enforce zero compatibility on i =  $i j = $j"
            mc[i, j][:, column] .*= 0
        end
    else
        # Set zero-valued non-dominant corners equal to dominant corner value
        @inbounds for j in axes(mc, 2), i in axes(mc, 1)
            i == dominant_i && j == dominant_j && continue
            if norm(mc[i, j][:, column]) < minnorm
                @debug "enforce zero compatibility - nonzero on i =  $i j = $j"
                mc[i, j][:, column] .= dom
            end
        end
    end
    mc
end


"""
    enforce_column_compatibility_with_dominant!(mc, dominant_i, dominant_j)

If necessary for valid interpolation, switch major-minor column order for non-dominant interpolation square corners.
"""
function enforce_column_compatibility_with_dominant!(mc, dominant_i, dominant_j)
    # We're checking column 1 only. We could be unlucky, then?
    dom = mc[dominant_i, dominant_j][:, 1]
    @inbounds for j in axes(mc, 2), i in axes(mc, 1)
        i == dominant_i && j == dominant_j && continue
        if is_close_to_perpendicular(mc[i, j][:, 1], dom)
            @debug "column swap: i = $i j = $j dominant_i = $(dominant_i) dominant_j = $(dominant_j) mc[i, j][:, 1] = $(mc[i, j][:, 1]) dom = $dom" maxlog = 120
            swap_columns!(mc, i, j)
        end
    end
    mc
end

"""
    enforce_direction_compatibility_with_dominant!(mc, dominant_i, dominant_j, column)

Two almost equal bidirectional vectors can point in opposite directions.
Here, we simply drop the interpolation in those cases. Flipping sign at non-
dominant corners of the interpolation square might be more accurate and ought be
tested out.
"""
function enforce_direction_compatibility_with_dominant!(mc, dominant_i, dominant_j, column)
    dom = mc[dominant_i, dominant_j][:, column]
    @inbounds for j in axes(mc, 2), i in axes(mc, 1)
        i == dominant_i && j == dominant_j && continue
        if is_close_to_opposite(mc[i, j][:, column], dom)
            @debug "enforce direction compat: i = $i j = $j column = $column mc[i, j][:, column] = $(mc[i, j][:, column]) dom = $dom" maxlog = 120
            mc[i, j][:, column] .= dom
        end
    end
    mc
end


"""
    swap_columns!(mc, i, j)
"""
function swap_columns!(mc, i, j)
    # This isn't called very often
    @debug "Swap columns i = $i j = $j"
    firstc = mc[i, j][:, 1]
    mc[i, j][:, 1] .= mc[i, j][:, 2]
    mc[i, j][:, 2] .= firstc
    mc
end

"""
    is_close_to_perpendicular(v1::T, v2::T) where T <: MVector{2, Float64}
    is_close_to_perpendicular(dotprod)

60°–120° returns true. If v1 or v2 are zero, returns false.

This interprets the two first elements as a 2d vector. Ignores other elements.
"""
@inline function is_close_to_perpendicular(v1::T, v2::T) where T <: MVector{2, Float64}
    l1 = norm(v1)
    #@show l1
    l1 == 0 && return false
    l2 = norm(v2)
    #@show l2
    l2 == 0 && return false
    d = dot_product(v1 / l1, v2 / l1)
    #@show d
    isit = is_close_to_perpendicular(d)
    if ! isit && abs(d) < 0.8
        @debug "Not perpendicular, but abs(d) = $(abs(d))"
    end
    #throw("ok..")
    isit
end
@inline is_close_to_perpendicular(dotprod::Float64) = abs(dotprod) < 0.5


"""
    is_close_to_opposite(dotprod)

Direction change is 138° < direction change <  228°
"""    
@inline function is_close_to_opposite(v1::T, v2::T) where T <: MVector{2, Float64}
    l1 = norm(v1)
    l1 == 0 && return false
    l2 = norm(v2)
    l2 == 0 && return false
    d = dot_product(v1 / l1, v2 / l1)
    is_close_to_opposite(d)
end
@inline is_close_to_opposite(dotprod::Float64) = dotprod < -0.668


@inline function dot_product(v1::T, v2::T) where T <: MVector{2, Float64}
    v1[1] * v2[1] + v1[2] * v2[2]
end



#########
# Various
######### 

"""
    normalize_or_zero!(K::TENSORMAP)
    normalize_or_zero!(v) where T<: SizedVector
    normalize_or_zero!(K::TENSORMAP, minnorm)
    normalize_or_zero!(v::T, minnorm)

Alternative to `normalize!` for those times when you want to make exceptions for zero.
The restrictive typing here errors on non-mutating calls like

```
julia> K
2×2 MMatrix{2, 2, Float64, 4} with indices SOneTo(2)×SOneTo(2):
 -1.0           1.12188e-14
 -5.59488e-12  -0.00200513

julia> normalize_or_zero!(K[:, 1])        # BAD CALL
ERROR: MethodError: no method matching normalize_or_zero!(::MVector{2, Float64})

julia> normalize_or_zero!(view(K, :, 1)); # MUTATING CALL

julia> K[:,1]
2-element MVector{2, Float64} with indices SOneTo(2):
 -1.0
 -5.594877687500939e-12
```
"""
function normalize_or_zero!(K::TENSORMAP, minnorm)
    normalize_or_zero!(view(K, :, 1), minnorm)
    normalize_or_zero!(view(K, :, 2), minnorm)
    K
end
function normalize_or_zero!(K::TENSORMAP)
    normalize_or_zero!(view(K, :, 1))
    normalize_or_zero!(view(K, :, 2))
    K
end
normalize_or_zero!(v::T) where T<: SizedVector =  normalize_or_zero!(v, MAG_EPS)
function normalize_or_zero!(v::T, minnorm) where T<: SizedVector
    mag = norm(v)
    if mag < minnorm
        v .= 0.0
    else
        v ./= mag
    end
    v
end



function enforce_minnorm!(K::TENSORMAP, minnorm::Float64)
    enforce_minnorm!(view(K, :, 1), minnorm)
    enforce_minnorm!(view(K, :, 2), minnorm)
    K
end
function enforce_minnorm!(v::T, minnorm) where T<: SizedVector
    if norm(v) < minnorm
        v .= 0.0
    end
    v 
end

#########################
# Derived properties of K
#########################

"""
    signed_curvature_values(K)
    --> {Float64, Float64}
"""
function signed_curvature_values(K)
    s1 = norm(K[:, 1]) * (is_bidirec_vect_positive(K[:, 1]) ? 1 : -1)
    s2 = norm(K[:, 2]) * (is_bidirec_vect_positive(K[:, 2]) ? 1 : -1)
    s1, s2
end

function rank_123_ranks(c1, c2)
    (c1 - 1) * (8 - c1) ÷ 2 + (c2 - c1 +1)
end


"""
    convexity_rank(s::Float64, flatval)
    --> Int 1-3
    convexity_rank(K, flatval)
    --> Int 1-6

Lexicographic rank, not "matrix rank".

# Arguments

- `s`       : signed curvature component along principal direction
- `flatval` : Flat is when  `abs(curvature) < flatval`
- `K` is a 2 x 2 matrix, format defined by `principal_curvature_components!`


| Signed curvature component | Rank   | Description |
| -------------------------- | -------| ----------- |
|            s > flatval     | 1      |  Convex     |
| -flatval < s < flatval     | 2      |  Flat       |
|            s < - flatval   | 3      |  Concave    |


Rank of principal curvature components (`K`):


| Major rank | Minor rank | Rank   | Description       |
| ---------- | ---------- | ------ | ----------------  |
|  1         | 1          |  1     | Convex  - convex  |
|  1         | 2          |  2     | Convex  - flat    |
|  1         | 3          |  3     | Convex  - concave |
|  2         | 2          |  4     | Flat    - flat    |
|  2         | 3          |  5     | Flat    - Concave |
|  3         | 3          |  6     | Concave - concave |

"""
convexity_rank(s::Float64, flatval) = s > flatval ? 1 : s < -flatval ? 3 : 2
function convexity_rank(K, flatval)
    # K is ordered: Coloumn 1 is major, column 2 is minor.
    s1, s2 = signed_curvature_values(K)
    @assert s1 >= s2 "s1 $s1 >= s2 $s2"
    c1 = convexity_rank(s1, flatval)
    c2 = convexity_rank(s2, flatval)
    # Ordered combination no., 1-6
    rank_123_ranks(c1, c2)
end