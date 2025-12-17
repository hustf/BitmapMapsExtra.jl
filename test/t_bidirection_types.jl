using Test
using BitmapMapsExtras
using BitmapMapsExtras.TestMatrices
using BitmapMapsExtras: BidirectionOnGrid, BidirectionInDomain, BidirectionAtXY, Domain
using BitmapMapsExtras: AbstractXYFunctor, SelectedVec2AtXY
using BitmapMapsExtras: 𝐊!, norm, reset!, 𝐊ᵤ!, mark_at!, RGB


#################
# Low-level tests
#################

@testset "Constant curvature" begin
    radius = TestMatrices.r
    z = z_cylinder(0.0)[300:500, 100:300];
    @testset "BidirectionOnGrid" begin
        bog = BidirectionOnGrid(𝐊!, z)
        pt = (100, 100)
        @test isapprox(bog(pt...), [0 0; 0   -1 / radius],  atol = 2e-4)
        @inferred bog(pt...)
    end
    @testset "BidirectionInDomain" begin
        bid = BidirectionInDomain(𝐊!, z; minnorm = 1e-4, normalize = false)
        # Exact on grid
        @test isapprox(bid(100.0, 100.0), [0 0; 0 -1 / radius],  atol = 2e-6)
        # 
        @inferred bid(100.0, 100.0)
        # Interpolated
        @test isapprox(bid(100.2, 100.0), [0 0; 0 -1/radius],  atol = 2e-6)
        @test isapprox(bid(100.7, 100.0), [0 0; 0 -1/radius],  atol = 2e-6)
        @test isapprox(bid(100.0, 100.2), [0 0; 0 -1/radius],  atol = 2e-6)
        @test isapprox(bid(100.0, 100.7), [0 0; 0 -1/radius],  atol = 2e-6)
        @test isapprox(bid(100.7, 100.7), [0 0; 0 -1/radius],  atol = 2e-6)
        @test isapprox(bid(100.7, 195.7), [0 0; 0 -1/radius],  atol = 2e-6)
    end
end
@testset "BidirectionAtXY" begin
    # Smoothly varying curvature, largest at centre lines
    z = z_paraboloid()
    baxy = BidirectionAtXY(𝐊!, z, true)
    @test baxy(3.0, 3.0) ≈ [0.0003794227888053568, 0.0008893924914720936]
    vbu = [norm(baxy(100.0, y)[:, 1]) for y in 50:0.2:53]
    @test sort(vbu) == vbu
    @test length(unique(vbu)) == length(vbu)
    vbu = [norm(baxy(x, 100.0)[:, 1]) for x in 50:0.2:53]
    @test sort(vbu, rev = true) == vbu
    @test length(unique(vbu)) == length(vbu)
    # Outside domain: zero everywhere, but no error thrown.
    # This is because we find it hard to restrict DiffEq solvers
    # from probing solutions outside the domain.
    @test baxy(2.99, 3.0) ≈ [0, 0]
    @test baxy(997.0, 997.0) ≈ [0.0003794227888056931, 0.0008893924914691132]
    @test baxy(997.01, 997.0) ≈ [0,  0]
    @test baxy(997.00, 997.01) ≈ [0,  0]
    @test baxy(892.0, 960.0) ≈  [0.0002916029305382096, 0.0009477262727801998]
    baxy.major[] = ! baxy.major[]
    @test baxy(892.0, 960.0) ≈ [-0.0012225868641178692, -0.0002940188601622608]
end

#@testset "Length of interpolated vector" begin
    # Smoothly varying curvature, largest at centre lines
    baxy = BidirectionAtXY(𝐊!, z_paraboloid(), true)
    # Domain border check.
    @test baxy(2.0, 2.0) == [0.0, 0.0]
    @test norm(baxy(3.0, 3.0)) ≈ 0.0009669439779861842
    vl = [norm(baxy(y, y)) for y in 2.0:0.2:3.0]
    @test all(iszero.(vl[1:5]))
    @test vl[6] ≈ 0.0009669439779861842
    # Interpolation of vectors in diverging region preserve length
    baxy = BidirectionAtXY(𝐊!, z_wavy(), true)
    img = background(baxy)

    rngx = 350:0.2:450
    vxy = [(x, -400 + 3x) for x in rngx]
    vpts = map(vxy) do (x,y)
        j = Int(round(x))
        i = baxy.negy(Int(round(y)))
        CartesianIndex((i,j))
    end
    mark_at!(img, vpts, 5, "in_circle")

    vl = [norm(baxy(xy...)) for xy in vxy]
    sort(vl, rev = true) == vl
#end

#=
# Dev speed
using BenchmarkTools
function foo(baxy)
    s = 0.0
    for x in 3:0.2:55, y in 3:0.2:55
        s += norm(baxy(x, y))
    end
    s
end
baxy = BidirectionAtXY(𝐊ᵤ!, z_paraboloid(), true)
@test foo(baxy) ≈ 68120.99400731127
@test baxy.K == [0.35340759825711493 -0.9690324166527553; 0.9354694380331928 -0.2469335446554816]
#0.237623 seconds (462.80 k allocations: 14.124 MiB)
@time foo(baxy) # maybe in precompilation
# 238.816 ms (462803 allocations: 14.12 MiB)
@btime foo(baxy)
=#




@testset "SelectedVec2AtXY" begin
    # Smoothly varying curvature, largest at centre lines
    z = z_paraboloid()
    saxy = SelectedVec2AtXY(𝐊!, z, true, false)
    @test saxy isa AbstractXYFunctor
    @test saxy(3.0, 3.0) ≈ [0.0003794227888056931, 0.0008893924914691132]
    vu = [norm(saxy(100.0, y)[:, 1]) for y in 50:0.2:53]
    @test sort(vu) == vu
    @test length(unique(vu)) == length(vu)
    vu = [norm(saxy(x, 100.0)[:, 1]) for x in 50:0.2:53]
    @test sort(vu, rev = true) == vu
    @test length(unique(vu)) == length(vu)
    @test size(saxy) == size(z)
    # Outside domain: zero everywhere, but no error thrown.
    # This is because we find it hard to restrict DiffEq solvers
    # from probing solutions outside the domain.
    @test saxy(2.99, 3.0) ≈ [0, 0]
    @test saxy(997.0, 997.0) ≈ [0.0003794227888056931, 0.0008893924914691132]
    @test saxy(997.01, 997.0) ≈ [0, 0]
    @test saxy(997.00, 997.01) ≈ [0, 0]
    @test saxy(892.0, 960.0) ≈ [0.0002916029305382096, 0.0009477262727801998]
    saxy.flip[] = ! saxy.flip[]
    @test saxy(892.0, 960.0) ≈ [-0.0002916029305382096, -0.0009477262727801998]
    saxy.baxy.major
    saxy.flip
    @test saxy.flip[]
    reset!(saxy)
    @test ! saxy.flip[]
end
