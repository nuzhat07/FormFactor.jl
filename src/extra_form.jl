
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

function form(n1::Int64,l1::Int64,m1::Int64,n2::Int64,l2::Int64,ρ::Float64,q::Float64)
	if n1>l1>=m1 && n2>l2>=m1
		Sum = 0
		for l in abs(l1-l2):2:l1+l2
			for k in 0:n2-l2-1
				for k1 in 0:n1-l1-1
					Sum += fun(n1,l1,m1,n2,l2,k,k1,l,ρ,q)
				end
			end
		end
		return Sum
	else
		"Incorrect Combination"
	end
end
