module Gaussian

using StaticArrays

export PrimitiveCartesianGaussian

struct PrimitiveCartesianGaussian
    ζ::Float64
    R::SVector{3, Float64}
    cartesianExponents::SVector{3, Int}
end

end
