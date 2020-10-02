module BasisFunctions

using Common, StaticArrays, OffsetArrays

export Gaussian, cartesianToSphericalMxs

struct Gaussian
    ζ::Float
    R::Coordinate
end

function positiveDiagonal(l::Int8, m::Int8)
    δl0 = l == 0 ? 1 : 0
    return sqrt(2 ^ δl0 * (2*l + 1) / (2 * l + 2)) 
end

function cartesianToSphericalMxs(lmax::Int)
    nang = (lmax+1)^2
    ncart = (lmax+1)*(lmax+2)*(lmax+3)÷6
    M = OffsetArray(Array([Matrix(zeros(Float64, ((l+1)*(l+2)÷2, 2*l+1))) for l=0:lmax]), -1)
    #  M = zeros(Float, (ncart, nang))
    for l=0:lmax #  Over common angular momentum quantum number
        iang = 0
        println()
        for m=-l:l  #  Spherical quantum number m
            iang += 1
            vm = m < 0 ? 0.5 : 0
            am = abs(m)
            vmax = Int(floor(am/2-vm))
            mfact = m == 0 ? 0.5 : 1
            N = 1 / (2^am * factorial(l)) * sqrt(2 * mfact * (factorial(l+am) * factorial(l-am)))
            #  icart = l ^ 2 #  Take into account the cartesians of all the previous ls
            icart = 0
            for k=0:l #  Cartesian quantum number k
                for j=0:l-k #  Cartesian quantum number j
                    i=l-k-j #  Cartesian quantum number i
                    #  if m == 0
                        #  println("x,y,z:",i,j,k)
                    #  end
                    icart += 1
                    for t=0:(l-am)÷2
                        kk = l - 2 * t - am
                        if kk ≠ k continue end
                        for u=0:t
                            for vv=0:vmax
                                v=vv+vm
                                ii = 2 * t + am - 2 * (u + v)
                                jj = 2 * (u + v)
                                if ii == i && jj == j && kk == k
                                    M[l][icart,iang] += (-1)^(t+vv) * 0.25 ^ t * N * binomial(l,t) * binomial(l-t,am+t) * binomial(t,u) * binomial(am,Int(2*v))
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return M
end

end
