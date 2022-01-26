"""
function fun(n1::Int64,l1::Int64,m1::Int64,a1::Float64,n2::Int64,l2::Int64,a2::Float64,q::Float64,k::Int64,k1::Int64,l::Int64)
	n = n1-l1-1
	δ = (n2 + a1*n1/a2)/(n1*n2)
	γ = l1+l2+k+3/2
	α = 2*l1+1
	ν = l+1/2
	N1 = 2^(l1+l2+k+k1-l+1)*(-1)^(m1+k+k1)*a1^(l2-l1+k)/(a2^(3/2+l2+k)*n1^(2+l1+k1)*n2^(2+l2+k))
	I1 = (1im*q)^l*(2*l+1)*gamma(ν+γ+k1+1)/(gamma(k+1)*gamma(n2-l2-k)*gamma(2*l2+k+2)*gamma(k1+1)*gamma(n1-l1-k1)*gamma(α+k1+1)*gamma(ν+1)*δ^(ν+γ+k1+1))
	A1 = √(π*gamma(n1-l1)*gamma(n2-l2)*gamma(n2+l2+1)*gamma(n1+l1+1)*(2*l1+1)*(2*l2+1))*wigner3j(l1, l2, l, 0, 0, 0)*wigner3j(l1, l2, l, m1, -m1, 0)*_₂F₁((ν+γ+k1+1)/2, (ν+γ+k1+2)/2, 1+ν, -q^2/δ^2)
	return N1*I1*A1
end
"""

 # for this special case m = m1-m2=0.



function fun(n1::Int64,l1::Int64,m1::Int64,a1::Float64,n2::Int64,l2::Int64,a2::Float64,q::Float64,k::Int64,k1::Int64,l::Int64)
	n = n1-l1-1
	δ = (a2*n2/a1+a1*n1)/(n1*a2*n2)
	γ = l1+l2+k+3/2
	β = 2/(a1*n1)
	α = 2*l1+1
	ν = l+1/2
	aa = (ν+γ+k1+1)/2
	bb = (ν+γ+k1+2)/2
	cc = 1+ν
	dd = -q^2/δ^2
	N1 = 4*π*4/(a2^(3/2)*(n1*n2)^2)*√(gamma(n1-l1)*gamma(n2-l2)/(gamma(n1+l1+1)*gamma(n2+l2+1)))*(2/(a1*n1))^l1*(2/(a2*n2/a1))^l2*√(π/(2*q))
	I1 = gamma(n2+l2+1)*(-2)^k*(-β)^k1*q^ν*gamma(n+α+1)*gamma(ν+γ+k1+1)/(gamma(k+1)*gamma(n2-l2-k)*gamma(2*l2+k+2)*(a2*n2/a1)^k*gamma(k1+1)*gamma(n-k1+1)*gamma(α+k1+1)*2^ν*gamma(ν+1)*δ^(ν+γ+k1+1))*_₂F₁(aa,bb,cc,dd)
	A1 = (1im)^l*√((2*l+1)*(2*l1+1)*(2*l2+1)*(2*l+1)/(4π*4π))*wigner3j(l1, l2, l, 0, 0, 0)*wigner3j(l1, l2, l, m1, -m1, 0)
	return N1*I1*A1
end

function form(n1::Int64,l1::Int64,m1::Int64,a1::Float64,n2::Int64,l2::Int64,a2::Float64,q::Float64)
	if n1>l1>=m1 && n2>l2>=m1
		Sum = 0
		for l in abs(l1-l2):2:l1+l2
			for k1 in 0:n2-l2-1
				for k in 0:n1-l1-1
					Sum += fun(n1,l1,m1,a1,n2,l2,a2,q,k,k1,l)
				end
			end
		end
		return Sum
	else
		"Incorrect Combination"
	end
end
