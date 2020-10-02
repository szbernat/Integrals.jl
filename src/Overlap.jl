module Overlap

using LinearAlgebra: dot
using StaticArrays: SVector
using OffsetArrays: OffsetArray
using Common: Float, Coordinate
using BasisFunctions: Gaussian

export buildCartesianOverlaps

"""
    gaussianOverlap(Rμ::Coordinate, ζμ::Float, Rν::Coordinate, ζν::Float)

Calculates the overlap integral ∫ μ⋅ν dr between the gaussian functions μ and ν having centers Rμ, Rν and exponents ζμ, ζν.
"""
function gaussianOverlap(Rμ::SVector{3, Float}, ζμ::Float, Rν::SVector{3, Float}, ζν::Float)
    dR = Rμ - Rν
    p = ζμ + ζν
    return sqrt(π / p) * exp(- ζμ * ζν / p * dot(dR, dR))
end

"""
    buildCartesianOverlaps(Rμ::SVector{3}, ζμ::Float, Rν::SVector{3}, ζν::Float, lμ::Int, lν::Int)

Calculates all of the primitive cartesian overlap integrals up to angular momentums lμ, lν between gaussians at centers Rμ, Rν with exponents ζμ, ζν.
"""
function buildCartesianOverlaps(Rμ::SVector{3, Float}, ζμ::Float, lμ::Int, Rν::SVector{3, Float}, ζν::Float, lν::Int)::AbstractArray
    p = ζμ + ζν
    twop = 0.5 / p
    S = OffsetArray(Array(zeros(Float, lμ+1, lν+1, 3)), -1, -1, 0)
    S[0,0,:] .= gaussianOverlap(Rμ, ζμ, Rν, ζν)
    for n=1:3 #  Over cartesian directions
        P = (ζμ * Rμ[n] + ζν * Rν[n]) / p
        XPμ = P - Rμ[n]
        XPν = P - Rν[n]
        if lμ > 0
            S[1,0,n] = XPμ * S[0,0,n]
        end
        if lν > 0
            S[0,1,n] = XPν * S[0,0,n]
        end
        if lμ > 0 && lν > 0
            S[1,1,n] = XPμ * S[0,1,n] + twop * S[0,0,n]
        end
        for i=1:lμ-1
            S[i+1,0,n] = XPμ * S[i,0,n] + twop * i * S[i-1,0,n]
        end
        for j=1:lν-1
            S[0,j+1,n] = XPν * S[0,j,n] + twop * j * S[0,j-1,n]
        end
        if lμ > 0
            for j=1:lν-1
                S[1,j+1,n] = XPν * S[1,j,n] + twop * (S[0,j,n] + j * S[1,j-1,n])
            end
        end
        for i =1:lμ-1
            for j=1:lν
                S[i+1,j,n] = XPμ * S[i,j,n] + twop * (i * S[i-1,j,n] + j * S[i,j-1,n])
            end
        end
    end
    return S
end

end
