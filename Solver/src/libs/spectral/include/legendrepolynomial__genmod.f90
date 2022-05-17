        !COMPILER-GENERATED INTERFACE MODULE: Tue May 17 10:04:26 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LEGENDREPOLYNOMIAL__genmod
          INTERFACE 
            FUNCTION LEGENDREPOLYNOMIAL(N,X) RESULT(L_N)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: X
              REAL(KIND=8) :: L_N
            END FUNCTION LEGENDREPOLYNOMIAL
          END INTERFACE 
        END MODULE LEGENDREPOLYNOMIAL__genmod
