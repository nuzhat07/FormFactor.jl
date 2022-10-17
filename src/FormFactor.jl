module FormFactor

using SpecialFunctions
using HypergeometricFunctions
using WignerFamilies
using FLoops
using AssociatedLegendrePolynomials



include("extra_form.jl")
include("independent_q.jl")

export mform,mfun,mfunq,mformq

end
