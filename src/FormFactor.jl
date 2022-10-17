module FormFactor

using SpecialFunctions
using HypergeometricFunctions
using WignerFamilies
using FLoops
using AssociatedLegendrePolynomials



include("extra_form.jl")
include("independent_q.jl")

export mform,mfun,mfunq,mformq
mformq(10,9,1,9,1,1,1.0,0.66,π/2)
end
