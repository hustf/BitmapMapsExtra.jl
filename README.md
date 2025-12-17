# BitmapMapsExtra

Experiments on top of the BitmapMaps pipeline. 

## Graphic output

BitmapMapsExtras can, based on an elevation matrix:

- Draw streamlines. Most typically downhill streamlines, but also
  lines of curvature. Individual lines are solutions to differential equations.
- Draw colorful and transparent glyphs:
  * 2d vector
  * 180° symmetrical 2d vector
  * Frames of reference (Darboux), local planes
- Place glyphs tightly packed respecting size, on grids, or random distributions
- Paint regions by convexity: 
  * Convex - convex
  * Convex - flat
  * Convex - concave
  * Flat - flat
  * Flat - concave
  * Concave - concave

## Differential geometry functions

- Normal vector projection, 𝐧ₚ
- Principal curvatures components, 𝐊

Both sample the elvation matrix, are allocation-free and fast. They have normalized-length versions, which is useful for drawing streamline plots. Both sample from a window of 5x5 in an elevation matrix and are defined on
discrete pixels.

Similar concepts, say stress, strain or creep tensors, could be plugged in to make similar graphics.


## `AbstractXYFunctor` and `AbstractIJFunctor`

Functors combine a function and its input matrix. The concrete types are

- `Vec2AtXY` and `Vec2OnGrid`
- `BidirectionAtXY` and `BidirectionOnGrid` 
- `SelectedVec2AtXY`


## `AbstractGlyphSpec` and `AbstractCurveSpec`

...for specifying visual appearance


## Some loose notes

`Vec2AtXY` enables integration along lines of descent or ascent. Visually, that
traces out a streamline. If we are only interested in the streamline shape for graphical purposes, we can integrate `𝐧ₚᵤ` (the direction of descent stripped from its steepness value) instead of `𝐧ₚ`. The result of integration along such a path is simple the elevation difference between start and finish. The direction of `𝐧ₚ` is always downhill, but we can move uphill by reversing the direction of 'time'.

Lines of curvature are made by picking one four directions, moving in that direction, and then picking direction again. If we do that consistently, we trace out a streamline.
The choice is between lines of maximal or minimal curvature, and in which of two directions on each of those.

 The result of integration along a line of curvature is change in steepness angle. However, since we're dealing with 2.5D height maps, this restricts the range to <-π, π>. By intuition, any streamline ought to be cyclic, or end in a point where we can't consistently move forward. For example, that might be a planar region.

Both types of streamline, they are more interesting in groups. For example, downhill streamlines collect in valleys, uphill streamlines collect in ridges and local summits. 

A collection of streamlines starting points is an ensemble. An ensemble might be a grid, uniformly random, or based on e.g. divergence (as in fluid flow) or previously drawn [hachures](https://andywoodruff.com/blog/hachures-and-sketchy-relief-maps/) in an iterative approach. 

We're using [OrdinaryDiffEq](https://docs.sciml.ai/OrdinaryDiffEq/stable/) to make these streamlines. Since our functors, or objects, are memory allocation free, we can find streamlines very effectively, even without using computational parallelism.

## Current progress

Version 0.1.5
Added paint_convexity and paint_steepness. Some more test geometry. 
Working on improvements to curvature lines, but real consistency still
eludes us. Tests in flux.

Version 0.1.4

Refactor functions in files.

Version 0.1.3

Improve painting by convexity type.

Version 0.1.2

Move graphics code to `DrawAndSpray.jl`
Drop troublesome dependency on `ColorBlendModes`.
Drop dependency on `BitmapMaps.jl`, temporarily in tests.

Version 0.0.15

First joint and working interface for streamlines:

- `plot_streamlines!(img, fxy::AbstractXYFunctor, pts; stroke = Stroke(), odekws...)`

The "functor" type names have been modified to highlight `vector` vs `bidirection`.

- Vec2AtXY{F, T}               <: AbstractXYFunctor
- SelectedVec2AtXY{F, T, LC}   <: AbstractXYFunctor
- Vec2OnGrid{F, T}             <: AbstractIJFunctor (for glyphs)
- BidirectionOnGrid{F, T, LC}  <: AbstractIJFunctor (for glyphs)

The streamlines for `SelectedVec2AtXY` has some ugly internals, but are consistently working. Revisions expected here.

In tests: `is_hash_stored` simplifies updating. All tests passing.

Version 0.0.14

Fully functional public API for all glyphs:
- `plot_glyphs!`
- `pack_glyphs!` 
- `indices_on_grid`

Glyph packing passes all tests.

- `plot_streamlines!` and `plot_bidirec_streamlines!` might also benefit from some similar, e.g. AbstractStreamlineSpec. Tests are currently broken, as we are dropping 
the temporary interface `grid_fcall` in 'common.jl' and will be working on a new streamlines interface.

Version 0.0.13

Added more crash tests for 'pack_glyphs'. Mostly working, but static dispatch on tensor glyph directions would be useful. Currently, packing 2d vectors will trigger an error. We plan to revise the GSTensor type.

Version 0.0.12

Add 'pack_glyphs' and restructure plotting function hierarchy. So far working for vectors, but more restructuring needed.

Version 0.0.11 Add missed file

Version 0.0.10 drops dependency ImageTransformations. Custom sub-pixel interpolations uses less memory allocation.

Version 0.0.9 outlines (currently in `test/t_calculate_and_draw_glyphs.jl`) a method to densely pack tensormap glyphs.

Version 0.0.8 adds a secondary color to the tensor basis glyph spec. Now, a two-color glyph can be drawn in one operation and faster.

Version 0.0.7 pre-allocates the linear solver for curvature directions, and makes other minor speed improvements.

We are revising and improving on a user interface, i.e. the high-level functions. The somewhat unsatisfactory interfaces reside in `test/common.jl`. Some test files are not updated to the latest version of the interface. One idea for a high-level interface is to pass types like glyph specification into more general functions. `grid_fcall_with_background` etc. is not a nice way to do this. 

`tangent_basis` has some minor allocations we should get rid of, but it may add more fields to some types.

Glyph stacking is currently implemented, but in an inefficient manner. A new algorithm for stacking glyps is planned, but a requirement for that algorithm is to know
the user's glyph specification in the stacking algorithm. Next up is to implement `glyphspec_types.jl` in `test/t_calculate_and_draw_glyphs`.

