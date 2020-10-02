using BasisFunctions
using Overlap
using OffsetArrays

@testset "Overlap.jl" begin
    Rμ = SVector{3}(randn(3))
    Rν = SVector{3}(randn(3))
    ζμ = abs(randn())
    ζν = abs(randn())
    @test all(buildCartesianOverlaps(Rμ, ζμ, 0, Rν, ζν, 0)[0,0,:] .≈ 1.6392143992458958)
    @test begin
        S = view(buildCartesianOverlaps(Rμ, ζμ, 2, Rν, ζν, 3), :, :, 2)
        Target = OffsetArray([1.6392143992 -0.0126705658 0.9068452213 -0.0210272948;
                              0.1026258527  0.9059540195 0.0427569217  1.5035741369;
                              0.9131723512  0.1064784544 1.5065815043  0.1534566213], -1, -1)
        S ≈ Target
    end
end
