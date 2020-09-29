using Gaussian
using Common, StaticArrays

@testset "Gaussian.jl" begin
    μ = PrimitiveCartesianGaussian(42, [1.,2.,3.], [0,1,0])

    @test μ.R ≈ [1,2,3]
    @test μ.cartesianExponents == [0,1,0]
    @test μ.ζ ≈ 42

    @test_throws TypeError PrimitiveCartesianGaussian(1,[0,0,0],[1])
    @test_throws DimensionMismatch PrimitiveCartesianGaussian(1,[0.,0.,0.],[1])
end
