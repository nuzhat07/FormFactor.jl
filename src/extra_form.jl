function mform(n1::Int64,l1::Int64,m1::Int64,a1::Float64,n2::Int64,l2::Int64,m2::Int64,a2::Float64,q::Float64)

    # The function wigner3j_f produces a vector of wigner symbols where each index corresponds to a value of l between
    # abs(l1-l2) and l1+l2 (and assumes m = -m1-m2) meaning that:
    #    collect(wigner3j_f(l1, l2, m1, m2)) ≈ [wigner3j(l1, l2, l3, m1, m2, -m1-m2) for l ∈ abs(l1-l2):l1+l2]

    wigner3j_l1_l2__0_0 = wigner3j_f(l1, l2, 0, 0)
    wigner3j_l1_l2__m1_negm2 = wigner3j_f(l1, l2, m1, -m2)

    δ = (a2*n2/a1+a1*n1)/(n1*a2*n2)
    β = 2/(a1*n1)
    α = 2*l1+1
    dd = -q^2/δ^2
    N1 = (4π*4/(a2^(3/2) * (n1 * n2)^2)
          * √(gamma(n1-l1) * gamma(n2-l2)
              / (gamma(n1+l1+1) * gamma(n2+l2+1)))
          * (2/(a1 * n1))^l1 * (2/(a2 * n2/a1))^l2 * √(π/(2q)))

    if n1 > l1 >= m1 && n2 > l2 >= m1
        # @floop from the package Floops.jl is able to parallelize a loop for us. If the loop contains a 'reduction'
        # (in our case, the sum), we just tell it which variable is being reduced with the @reduce macro and it takes
        # care of initializing the reduction variables and will make sure the multi-threading does not introduce a 'race
        # condition' which could cause the sum to be incorrect.
        @floop for l in abs(l1-l2):2:l1+l2
            for k in 0:n2-l2-1
                for k1 in 0:n1-l1-1
                    if (l1+l2) % 2 == 0
                        @reduce Sum += mfun(n1,l1,m1,a1,n2,l2,m2,a2,q,k,k1,l, δ, β, α, dd, N1,
                                            wigner3j_l1_l2__0_0, wigner3j_l1_l2__m1_negm2) 
                    else
                        @reduce Sum += im * mfun(n1,l1,m1,a1,n2,l2,m2,a2,q,k,k1,l,  δ, β, α, dd, N1,
                                                 wigner3j_l1_l2__0_0, wigner3j_l1_l2__m1_negm2)
                    end
                end
            end
        end
        return N1 * Sum
    else
        "Incorrect Combination"
    end
end

function mfun(n1::Int64,l1::Int64,m1::Int64,a1::Float64,n2::Int64,l2::Int64,m2::Int64,a2::Float64,q::Float64,k::Int64,k1::Int64,l::Int64,
              δ, β, α, dd, N1,
              wigner3j_l1_l2__0_0, wigner3j_l1_l2__m1_negm2)
    m=0
    n = n1-l1-1
    γ = l1+l2+k+3/2
    ν = l+1/2
    aa = (ν+γ+k1+1)/2
    bb = (ν+γ+k1+2)/2
    cc = 1+ν

    I1 = (gamma(n2+l2+1) * (-2)^k * (-β)^k1 * q^ν * gamma(n+α+1)*gamma(ν+γ+k1+1)
          / (gamma(k+1) * gamma(n2-l2-k) * gamma(2*l2+k+2) * (a2*n2/a1)^k * gamma(k1+1)
             * gamma(n-k1+1) * gamma(α+k1+1) * 2^ν * gamma(ν+1) * δ^(ν+γ+k1+1)) * _₂F₁(aa,bb,cc,dd))

    if (l1+l2) % 2 == 0
        A1 = (-1)^(l/2+m1) * √((2*l+1) * (2*l1+1) * (2*l2+1) * (2*l+1) / (4π*4π))
    else
        A1 = (-1)^((l-1)/2+m1) * √((2*l+1) * (2*l1+1) * (2*l2+1) * (2*l+1) / (4π*4π))
    end
    return I1 * A1 * wigner3j_l1_l2__0_0[l] * wigner3j_l1_l2__m1_negm2[l]
end
