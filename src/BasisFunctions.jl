module BasisFunctions

using Common, StaticArrays

export Gaussian

struct Gaussian
    ζ::Float
    R::Coordinate
end

function positiveDiagonal(l::Int8, m::Int8)
    δl0 = l == 0 ? 1 : 0
    return sqrt(2 ^ δl0 * (2*l + 1) / (2 * l + 2)) 
end

end
