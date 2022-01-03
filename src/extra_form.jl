kron(x, y) = ==(x, y)
function kron(x, y, z...)
	!kron(x, y) && return false
	for i in z
		!kron(x, i) && return false
	end
	return true
end

function form(n1::Int64,l1::Int64,m1::Int64,a1::Int64,n2::Int64,l2::Int64,m2::Int64,a2::Int64,q::Float64)
	m = 0 # for this special case m = m1-m2=0.
	n = Int(n1-l1-1)

	    function fun(k,k1,l)
	        γ = l1+l2+k1+3/2
	        δ = ((a2*n2)+(a1*n1))/(n1*n2*a2*a1)
	        β = 2/(a1*n1)
	        μ = q
	        α = 2*l1+1
	        ν = l+1/2
	        aa = (ν+γ+k+1)/2
	        bb = (ν+γ+k+2)/2
	        cc = 1+ν
	        dd = -(μ^2/δ^2)
	        N1 = 4*pi*(2/(a1^(3/2)*n1^2))*(2/((a2/a1)^(3/2)*n2^2))*√(factorial(big(n1-l1-1))*factorial(big(n2-l2-1))/(factorial(big(n1+l1))*factorial(big(n2+l2))))*(2/(a1*n1))^l1*(2*a1/(a2*n2))^l2*√(pi/(2*q))*kron(m,0)

	        Il = factorial(big(n2+l2))*(-1)^k1*2^k1/(factorial(big(k1))*factorial(big(n2-l2-1-k1))*factorial(big(2*l2+k1+1))*(a2)^k1*n2^k1)*((-β)^k*μ^ν*gamma(n+α+1)*gamma(ν+γ+k+1)/(factorial(big(k))*factorial(big(n-k))*gamma(α+k+1)*2^ν*gamma(ν+1)*δ^(ν+γ+k+1)))*_₂F₁(aa,bb, cc,dd)
	        Al = (complex(1))^l*√((2*l+1)/(4*pi))*(-1)^(m2+m)*√((2*l1+1)*(2*l2+1)*(2*l+1)/(4*pi))*wigner3j(l1, l2, l, 0, 0, 0)*wigner3j(l1, l2, l, m1, -m2, -m)
	        return N1*Il*Al
		end


		X = Float64[]
		Sum = 0

	    for p in range(abs(l1-l2),l1+l2)
	        for i in range(0,n)
	            for j in range(0,n2-l2-1)
	                Sum += fun(i,j,p)
					Sum =  convert(AbstractFloat, Sum)
				end
			end
		end
		push!(X,Sum)
		Y = Float64.(Tuple(X))
        @printf "%.5f" Y[1]
end
