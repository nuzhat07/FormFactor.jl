# FormFactor

[![Build Status](https://travis-ci.com/nuzhat07/FormFactor.jl.svg?branch=main)](https://travis-ci.com/nuzhat07/FormFactor.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/nuzhat07/FormFactor.jl?svg=true)](https://ci.appveyor.com/project/nuzhat07/FormFactor-jl)
[![Coverage](https://codecov.io/gh/nuzhat07/FormFactor.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nuzhat07/FormFactor.jl)
[![Coverage](https://coveralls.io/repos/github/nuzhat07/FormFactor.jl/badge.svg?branch=main)](https://coveralls.io/github/nuzhat07/FormFactor.jl?branch=main)

**Introduction**

This is a repository for atomic form factor calculation. The function "mform(n1::Int64,l1::Int64,m1::Int64,n2::Int64,l2::Int64,ρ::Float64,q::Float64)" will calculate the form factor and the input of this function are  
 * n1, l1, m1 are the quantum numbers of the initial particle and n2, l2 are the quantum numbers of the final particle.
 * ρ is the ratio of Bohr radii of the initial and final particle. We are using the length in the units of Bohr radius of initial particle.
 * q is the transferred momentum in atomic unit. Atomic units: qatomic units = q/ (μcα), where α is the fine structure constant, μ is the reduced mass of the hydrogen-like atom.


 **Installation**

 To install this package one has to do

 * In the Julia REPL package mode `pkg> add https://github.com/nuzhat07/FormFactor.jl.git`
 * `julia> import FormFactor`
 * `julia> using FormFactor`
 * `julia> mform(1,0,0,1,0,1.0,0.66)`
 
 **Fortran code**
 The author of [[1]](#1) implemented a different algorithm in Fortran to calculate atomic form factor. To compare the computation time, we have changed the *ffdiff.f* file of the package `https://elsevier.digitalcommonsdata.com/datasets/ddy58t53dz/1` and computed the computation time. 
 ```fortran
*
*
* *** Test program to calculate electron-hydrogen scattering cross
* *** sections.
*
************************************************************************
*
      IMPLICIT NONE
      INTEGER N, L, M
      DOUBLE PRECISION DIFOFA
*
      WRITE(6,*) '                                      '
      WRITE(6,*) '** ******************************** **'
      WRITE(6,*) '**                                  **'
      WRITE(6,*) '**  This Output File gives us the   **'
      WRITE(6,*) '**                                  **'
      WRITE(6,*) '**  absolute value of the result    **'
      WRITE(6,*) '**                                  **'
      WRITE(6,*) '**  of subtracting 1. to the form   **'
      WRITE(6,*) '**                                  **'
      WRITE(6,*) '**  form factor of a state at the   **'
      WRITE(6,*) '**                                  **'
      WRITE(6,*) '**  origin. The expected result is  **'
      WRITE(6,*) '**                                  **'
      WRITE(6,*) '**  0.                              **'
      WRITE(6,*) '** ******************************** **'
      WRITE(6,*) '                                      '
      WRITE(6,*) '--------------------------------------'
      WRITE(6,*) '|  N |  L |  M |  | FFactor(0.66)|   |'
      WRITE(6,*) '--------------------------------------'
      DO 1 N=1,30
            DO 2 L=0,N-1
              DO 3 M=0,L
               WRITE(6,101) N,L,M,ABS(DIFOFA(N,L,M,N,L,M,0.66D0))
          3            CONTINUE
          2          CONTINUE
          1       CONTINUE
      WRITE(6,*) '--------------------------------------'
*
 101  FORMAT(' | ',3(I2,' | '),ES11.2E3,'      | ')
      END
```



The author of [[2]](#2) implemented the same algorithm in Fortran for equal masses as we. To compare the computation time with this Fortran code, we have used the *crsmumu.f* file from the MuMuPy package `https://data.mendeley.com/datasets/nr6y34yg29/1` [[2]](#2) and written a *main.f* file.
 
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
## References
<a id="1">[1]</a> 
C. S. Ríos and J. S. Silva, 
An implementation of atomic form factors, 
Computer Physics Communications 151, (2003).

<a id="2">[2]</a> 
A. Uskov, A. Alizzi, and Z. Silagadze, 
MuMuPy: A dimuonium-matter interaction calculator, 
Computer Physics Communications 276 (2022). 
 
