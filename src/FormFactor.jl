module FormFactor

using SpecialFunctions
using HypergeometricFunctions
using WignerFamilies
using FLoops



include("extra_form.jl")

export mform,mfun
@time for n1 ∈ 1:50
    for l1 ∈ 0:n1-1
               for m1 ∈ 0:l1
                   print("$n1 " , "$l1 " , "$m1 " , mform(n1, l1, m1, n1, l1, 1.0, 0.66),"\n")
               end
           end
       end
mform(10, 8, 7, 8, 7, 1.0, 0.66)
@time mform(5,4,3,4,3,1.0,0.66)
end
