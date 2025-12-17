using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: principal_curvature_components, allocations_curvature
using BitmapMapsExtras: principal_curvature_components!, VΦ, sample_at_float
using BitmapMapsExtras: sampled_second_derivative_in_direction, norm
using BitmapMapsExtras: normalize_or_zero!
@testset "Curvature and normalization" begin
    Ω = CartesianIndices((-2:2, -2:2))
    M = z_cylinder(0)[TestMatrices.I0 + CartesianIndex(100, 100) .+ Ω] 
    CM = principal_curvature_components(z_cylinder(0), TestMatrices.I0)
    @test isapprox(CM, [-9.372522023933059e-5 0; 0 -0.0020051318513509423])
    @test isapprox(CM[1,1], 0, atol = 1e-4)
    @test isapprox(CM[2,2], -1 / TestMatrices.r, atol = 1e-5)
    # How not to pass arguments to the mutating function
    lv1 = norm(CM[:, 1])
    lv2 = norm(CM[:, 2])
    @test norm(normalize_or_zero!(CM[:, 1])) == 1
    @test norm(CM[:, 1]) == lv1
    # How to pass arguments to the mutating function
    normalize_or_zero!(view(CM, :, 1))
    @test norm(CM[:, 1]) == 1
end