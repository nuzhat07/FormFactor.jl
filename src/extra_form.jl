kron(x, y) = ==(x, y)
function kron(x, y, z...)
	!kron(x, y) && return false
	for i in z
		!kron(x, i) && return false
	end
	return true
end

function fun(n1,l1,m1,a1,n2,l2,m2,a2,q,k,k1,l)
	m = 0 # for this special case m = m1-m2=0.
	n = n1-l1-1
	γ = l1+l2+k+3/2
	δ = (a2*n2 + a1*n1)/(n1*n2*a2)
	β = 2/(n1)
	μ = q
	α = 2*l1+1
	ν = l+1/2
	aa = (ν+γ+k1+1)/2
	bb = (ν+γ+k1+2)/2
	cc = 1+ν
	dd = -(μ^2/δ^2)
	N1 = ((-1)^(m2+m)*2^(2+l1+l2))/((a2/a1)^(3/2+l2)*n1^(2+l1)*n2^(2+l2))*kron(m,0)*√(pi/(2*q))*√(((2*l1+1)*(2*l2+1)*factorial(n1+l1)*factorial(n2+l2)*factorial(n1-l1-1))/factorial(n2-l2-1))

	Il = ((-β)^k1*gamma(ν+γ+k1+1))/(factorial(k1)*gamma(n-k1+1)*gamma(α+k1+1)*δ^(ν+γ+k1+1))*_₂F₁(aa,bb,cc,dd)
	Al = ((complex(1))^l*(-1)^k*2^k*μ^ν*(2*l+1))/(2^ν*gamma(ν+1)*n2^k*(a2/a1)^k*factorial(k)*factorial(2*l2+k+1))*wigner3j(l1, l2, l, 0, 0, 0)*wigner3j(l1, l2, l, m1, -m2, -m)
	return N1*Il*Al
end

function form(n1::Int64,l1::Int64,m1::Int64,a1,n2::Int64,l2::Int64,m2::Int64,a2,q::Float64)


		#X = Float64[]
		Sum = 0

	    for p in range(abs(l1-l2),l1+l2)
	        for i in range(0,n2-l2-1)
	            for j in range(0,n1-l1-1)
	                Sum += fun(n1,l1,m1,a1,n2,l2,m2,a2,q,i,j,p)
					Sum =  convert(Float64, Sum)
				end
			end
		end
		#push!(X,Sum)
		#Y = Float64.(Tuple(X))
        return Sum #@printf "%.5f" Sum #Sum #Y[1] #@printf "%.5f" Y[1] #
end
