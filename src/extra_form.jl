function fun(n1::Int64,l1::Int64,m1::Int64,a1::Float64,n2::Int64,l2::Int64,m2::Int64,a2::Float64,q::Float64,k::Int64,k1::Int64,l::Int64)
	m=0
	n = n1-l1-1
	δ = (n2 + a1*n1/a2)/(n1*n2)
	γ = l1+l2+k+3/2
	β = 2/n1
	α = 2*l1+1
	ν = l+1/2
	N1 = 2^(l1+l2+k+k1-l+1)*(-1)^(m2+m+k+k1)/(a2^(3/2+l2)*n1^(2+l1+k1)*n2^(2+l2+k))*a1^(l2-l1)*(a1/a2)^k
	I1 = gamma(ν+γ+k1+1)/(gamma(k+1)*gamma(n2-l2-k)*gamma(2*l2+k+2)*gamma(k1+1)*gamma(n-k1+1)*gamma(α+k1+1)*gamma(ν+1)*δ^(ν+γ+k1+1))
	A1 = (1im*q)^l*(2*l+1)*√(π*gamma(n1-l1)*gamma(n2-l2)*gamma(n2+l2+1)*gamma(n1+l1+1)*(2*l1+1)*(2*l2+1))*wigner3j(l1, l2, l, 0, 0, 0)*wigner3j(l1, l2, l, m1, -m2, -m)*_₂F₁((ν+γ+k1+1)/2, (ν+γ+k1+2)/2, 1+ν, -q^2/δ^2)
	return N1*I1*A1
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
