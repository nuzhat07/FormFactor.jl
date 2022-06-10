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
 
 **Fortran code**
 
 We have used the *crsmumu.f* file from the MuMuPy package `https://data.mendeley.com/datasets/nr6y34yg29/1` and written a *main.f* file to compare the computation time.
 
 ```fortran
    PROGRAM FORM
      IMPLICIT NONE
      INTEGER n1,m1,l1,array1(2)
      REAL*8 Q, FTRANS
      LOGICAL loop
      EXTERNAL FTRANS
      PARAMETER (Q = 0.66D0)
      PARAMETER (array1 = (/0, 1/))
      WRITE(*,*)
      loop = .true.
      IF(loop.eqv..true.) THEN 
        DO n1=1,30
          DO l1=0,n1-1
            DO m1=0,l1
               WRITE(6,101) n1, l1, m1, ABS(FTRANS(array1,n1,l1,m1,n1,
     &       l1,m1,0.66D0))
            END DO
          END DO
       END DO
      ELSE
         WRITE(6,102) 2, 1, 1, 2, 1, 1, ABS(FTRANS(array1,2,1,
     &          1,2,1,1,0.66D0))    
      END IF
      IF(loop.eqv..true.) THEN
 101     FORMAT(' | ',3(I2,' | '),ES11.2E3,' | ')
      ELSE
 102     FORMAT(' | ',6(I2,' | '),ES11.2E3,' | ')
      END IF
      END PROGRAM FORM
```
 
 
