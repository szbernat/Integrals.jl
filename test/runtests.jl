using Test, Random

Random.seed!(271828)

include("Integrals_test.jl")
include("BasisFunctions_test.jl")
include("Overlap_test.jl")
