
"""
function fun(n1::Int64,l1::Int64,m1::Int64,n2::Int64,l2::Int64,k::Int64,k1::Int64,l::Int64,ρ::Float64,q::Float64)
	n = n1-l1-1
	γ = l1+l2+k+3/2
	δ = (n2 + n1*ρ)/(n1*n2)
	ν = l+1/2
	N1 = (-1)^(m1+k1+k) * 2^(k+k1+l1+l2-l+1) * ρ^(3/2+l2+k)*
    √(pi*(2*l1+1)*(2*l2+1)*factorial(n1+l1)*factorial(n2+l2)*factorial(n1-l1-1)/
        factorial(n2-l2-1)) / (n1^(2+l1)*n2^(2+l2))
	Il = gamma(ν+γ+k1+1)/(n1^k1*factorial(k1)*gamma(n-k1+1)*gamma(2*l1+k1+2)*δ^(ν+γ+k1+1))*
    _₂F₁((ν+γ+k1+1)/2, (ν+γ+k1+2)/2, 1+ν, -q^2/δ^2)
	Al = (1im*q)^l *(2*l+1)*wigner3j(l1, l2, l, 0, 0, 0)*
    wigner3j(l1, l2, l, m1, -m1, 0) / (n2^k*gamma(ν+1)*factorial(k)*factorial(2*l2+k+1))
	return N1*Il*Al
end


function fun(n1::Int64,l1::Int64,m1::Int64,a1::Float64,n2::Int64,l2::Int64,m2::Int64,a2::Float64,q::Float64,k::Int64,k1::Int64,l::Int64)
	m = 0 # for this special case m = m1-m2=0.
 	n = n1-l1-1
 	γ = l1+l2+k+3/2
 	δ = (a2*n2 + a1*n1)/(n1*n2*a2*a1)
 	#δ = (a2/a1*n2 + n1)/(n1*n2*a2/a1)
 	β = 2/n1
 	#μ = q
 	α = 2*l1+1
 	ν = l+1/2
 	aa = (ν+γ+k1+1)/2
 	bb = (ν+γ+k1+2)/2
 	cc = 1+ν
 	#dd = -q^2/δ^2
 	dd = -q^2/δ^2
 	N1 = ((-1)^(m2+m)*2^(2+l1+l2))/((a2/a1)^(3/2+l2)*n1^(2+l1)*n2^(2+l2))*√(pi/(2*q))*√(((2*l1+1)*(2*l2+1)*factorial(n1+l1)*factorial(n2+l2)*factorial(n1-l1-1))/factorial(n2-l2-1))
 	Il = ((-β)^k1*gamma(ν+γ+k1+1))/(factorial(k1)*gamma(n-k1+1)*gamma(α+k1+1)*δ^(ν+γ+k1+1))*_₂F₁(aa,bb,cc,dd)
 	#Al = (1im)^l*(-2)^k*q^ν*(2*l+1)/(2^ν*gamma(ν+1)*(n2*a2/a1)^k*factorial(k)*factorial(2*l2+k+1))*wigner3j(l1, l2, l, 0, 0, 0)*wigner3j(l1, l2, l, m1, -m2, -m)
 	Al = (1im)^l*(-2)^k*q^ν*(2*l+1)/(2^ν*gamma(ν+1)*(n2*a2/a1)^k*factorial(k)*factorial(2*l2+k+1))*wigner3j(l1, l2, l, 0, 0, 0)*wigner3j(l1, l2, l, m1, -m2, -m)
 	return N1*Il*Al
 end


function form(n1::Int64,l1::Int64,m1::Int64,a1::Float64,n2::Int64,l2::Int64,m2::Int64,a2::Float64,q::Float64)
	if n1>l1>=m1 && n2>l2>=m2
		Sum = 0
		for l in abs(l1-l2):l1+l2
			for k in 0:n2-l2-1
				for k1 in 0:n1-l1-1
					Sum += fun(n1,l1,m1,a1,n2,l2,m2,a2,q,k,k1,l)
				end
			end
		end
		return Sum
	else
		"Incorrect Combination"
	end
end
"""
function fun(n1::Int64,l1::Int64,m1::Int64,a1::Float64,n2::Int64,l2::Int64,m2::Int64,a2::Float64,q::Float64,k::Int64,k1::Int64,l::Int64)
	m=0
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
	A1 = (1im)^l*√((2*l+1)*(2*l1+1)*(2*l2+1)*(2*l+1)/(4π*4π))*wigner3j(l1, l2, l, 0, 0, 0)*wigner3j(l1, l2, l, m1, -m2, -m)
	return N1*I1*A1
end

function form(n1::Int64,l1::Int64,m1::Int64,a1::Float64,n2::Int64,l2::Int64,m2::Int64,a2::Float64,q::Float64)
	Sum = 0
	for l in abs(l1-l2):l1+l2
		for k in 0:n2-l2-1
			for k1 in 0:n1-l1-1
				Sum += fun(n1,l1,m1,a1,n2,l2,m2,a2,q,k,k1,l)
			end
		end
	end
	return Sum
end
