# FormFactor

[![Build Status](https://travis-ci.com/nuzhat07/FormFactor.jl.svg?branch=main)](https://travis-ci.com/nuzhat07/FormFactor.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/nuzhat07/FormFactor.jl?svg=true)](https://ci.appveyor.com/project/nuzhat07/FormFactor-jl)
[![Coverage](https://codecov.io/gh/nuzhat07/FormFactor.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nuzhat07/FormFactor.jl)
[![Coverage](https://coveralls.io/repos/github/nuzhat07/FormFactor.jl/badge.svg?branch=main)](https://coveralls.io/github/nuzhat07/FormFactor.jl?branch=main)

**Introduction**

This is a repository for atomic form factor calculation. The function "form" will calculate the form factor and the input of this function are  
 * n1, l1, m1 are the quantum numbers of the initial particle and n2, l2, m2 are the quantum numbers of the final particle.
 * a1 and a2 are the Bohr radius of the initial and final particle respectively. We are using the length in the units of Bohr radius of initial particle
 * q is the transferred momentum in atomic unit.


 **Installation**

 To install this package one has to do

 * In the Julia REPL package mode pkg> add 'github repository link of this package'
 * julia> import FormFactor
 * julia> using FormFactor
 * julia> form(1,0,0,1,1,0,0,1,0.66)
