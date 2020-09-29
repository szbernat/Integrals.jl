module Gaussian

using Common, StaticArrays

export PrimitiveCartesianGaussian

struct PrimitiveCartesianGaussian
    ζ::Float
    R::Coordinate
    cartesianExponents::SVector{3, Int}
end

end
