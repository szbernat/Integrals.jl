using BasisFunctions
using Overlap
using OffsetArrays

@testset "Overlap.jl" begin
    μ = Gaussian(randn(), randn(3))
    ν = Gaussian(randn(), randn(3))
    @test all(buildCartesianOverlaps(μ, ν, 0, 0)[0,0,:] .≈ 0.5050879913)
    @test begin
        S = view(buildCartesianOverlaps(μ, ν, 2, 3), :, :, 2)
        Target = OffsetArray([0.5050879913 -0.2868346941  0.3451125733 -0.4029503126;
                              0.2908210474  0.0170673400 -0.0082542143  0.1415094741;
                              0.3496716598  0.0112652498  0.1320693481 -0.0728283362], -1, -1)
        S ≈ Target
    end
end
