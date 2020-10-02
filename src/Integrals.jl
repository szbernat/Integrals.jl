module Integrals

using OffsetArrays: OffsetArray
using StaticArrays: SVector, MMatrix
using BasisFunctions: Gaussian, cartesianToSphericalMxs
using Common: Float
using Overlap: buildCartesianOverlaps

export OneElectronIntegral, OverlapIntegral, contractedIntegrals, assembleCartesianIntegrals

"""
    OneElectronIntegral

Abstract type that has all one electron integrals as its subtypes.
"""
abstract type OneElectronIntegral end

"""
    OverlapIntegral <: OneElectronIntegral

Abstract type to represent one electron overlap integrals.
"""
abstract type OverlapIntegral <: OneElectronIntegral end

"""
    contractedIntegrals(::Type{OverlapIntegral}, Rμ::SVector{3}, cμ::SVector, ζμ::SVector, Rν::SVector{3}, cν::SVector, ζν::SVector, lμ::Int, lν::Int)

Calculates the overlap integrals between contracted cartesian gaussians at centers Rμ, Rν with expansion coeffiencts cμ, cν and exponents ζμ, ζν, up to angular momentums lμ and lν.
"""
function contractedIntegrals(::Type{OverlapIntegral}, Rμ::SVector{3, Float}, cμ::SVector{Nμ, Float}, ζμ::SVector{Nμ, Float}, Rν::SVector{3, Float}, cν::SVector{Nν, Float}, ζν::SVector{Nν, Float}, lμ::Int, lν::Int) where {Nμ, Nν}
    contracted = OffsetArray(zeros(Float, lμ+1, lν+1, 3), -1, -1, 0)
    for i=1:Nμ, j=1:Nν
        contracted += cμ[i] * cν[j] * buildCartesianOverlaps(Rμ, ζμ[i], lμ, Rν, ζν[j], lν)
    end
    return contracted
end

"""
    assembleCartesianIntegrlas(::Type{<:Integral}, S::IntsInCartDirs, lμ, lν)

Assemble all (lμ+1)*(lμ+2)÷2 * (lν+1)*(lν+2)÷2 cartesian integrals from S which contains the integrals factorized into cartesian directions.
"""
function assembleCartesianIntegrals(::Type{OneElectronIntegral}, lμ::Int, lν::Int)
#  function assembleCartesianGTOIntegrals(::Type{OneElectronIntegral}, S::AbstractArray, lμ::Int, lν::Int)
    Nμ = (lμ+1)*(lμ+2)÷2
    Nν = (lν+1)*(lν+2)÷2
    cartesianGTO = MMatrix{Nμ, Nν, Float}(undef)
    nν=0
    S = contractedIntegrals(OverlapIntegral, SVector{3, Float}(randn(3)), SVector{10, Float}(randn(10)), SVector{10, Float}(abs.(randn(10))), SVector{3, Float}(randn(3)), SVector{6, Float}(randn(6)), SVector{6, Float}(abs.(randn(6))), lμ, lν)
    for kν=0:lν
        for jν=0:lν-kν
            iν=lν-kν-jν
            nν += 1
            nμ = 0
            for kμ=0:lμ
                for jμ=0:lμ-kμ
                    iμ=lμ-kμ-jμ
                    nμ += 1
                    cartesianGTO[nμ, nν] = S[iμ, iν, 1] * S[jμ, jν, 2] * S[kμ, kν, 3]
                end
            end
        end
    end
    return cartesianGTO
end

function overlap()
    l = 8
    Scart = assembleCartesianIntegrals(OneElectronIntegral, l, l)
    cartToSpher = cartesianToSphericalMxs(l)
    Sspher = transpose(cartToSpher[l]) * Scart * cartToSpher[l]
    return Sspher
end

end
