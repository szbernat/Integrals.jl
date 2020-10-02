module Common

using StaticArrays

export Float, Coordinate, Coordinate2, AngMomArray

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
    Coordinate(coords::NTuple{3, T}) where {T} = throw(TypeError(:Coordinate, "coords", Type{Float}, T))
end
Base.getindex(R::Coordinate, i::Int) = R.coords[i]
Tuple(R::Coordinate) = R.coords

#  struct AngMomArray{L, S, T}
    #  data::NTuple{S, T}
#  end

#  function AngMomArray{L, T}(data) where {L<:Tuple, T}
    #  S = (L.parameters[1]+1)^2
    #  AngMomArray{L, S, T}(data)
#  end
#  Base.getindex(a::AngMomArray{L, S, T}, l::Int, m::Int) where {L, S, T} = a.data[l^2+m+l+1]

end
