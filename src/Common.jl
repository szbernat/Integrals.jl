module Common

using StaticArrays

export Float, Coordinate, Coordinate2

Float = Float64

"""
    Coordinate(x::NTuple{3, Float64})
    Coordinate(x::AbstractArray{3, Float64})
    Coordinate(x1, x2, x3)

Represents the cartesian coordinates of a point in 3D space.
"""
struct Coordinate <: StaticArray{Tuple{3}, Float, 1}
    coords::NTuple{3, Float}
    Coordinate(coords::NTuple{3, Float}) = new(coords)
end
Base.getindex(R::Coordinate, i::Int) = R.coords[i]
Tuple(R::Coordinate) = R.coords

end
