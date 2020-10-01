using BasisFunctions
using Common, StaticArrays

@testset "BasisFunctions.jl" begin
    μ = Gaussian(42, [1.,2.,3.])

    @test μ.R ≈ [1,2,3]
    @test μ.ζ ≈ 42

    @test_throws TypeError Gaussian(1,[0,0,0])
    @test_throws DimensionMismatch Gaussian(1,[0.,0.])
end
