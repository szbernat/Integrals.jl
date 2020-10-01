module Overlap

using LinearAlgebra: dot
using OffsetArrays: OffsetArray
using Common: Float
using BasisFunctions: Gaussian

export buildCartesianOverlaps

function gaussianOverlap(μ::Gaussian, ν::Gaussian)
    dR = μ.R - ν.R
    p = μ.ζ + ν.ζ
    return sqrt(π / p) * exp(- μ.ζ * ν.ζ / p * dot(dR, dR))
end

function buildCartesianOverlaps(μ::Gaussian, ν::Gaussian, lμ::Int, lν::Int)::AbstractArray
    p = μ.ζ + ν.ζ
    twop = 0.5 / p
    S = OffsetArray(Array(zeros(Float, lμ+1, lν+1, 3)), -1, -1, 0)
    S[0,0,:] .= gaussianOverlap(μ, ν)
    for n=1:3 #  Over cartesian directions
        P = (μ.ζ * μ.R[n] + ν.ζ * ν.R[n]) / p
        XPμ = P - μ.R[n]
        XPν = P - ν.R[n]
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
