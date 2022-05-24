        !COMPILER-GENERATED INTERFACE MODULE: Tue May 24 12:21:47 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE USERDEFINEDSTATE1__genmod
          INTERFACE 
            SUBROUTINE USERDEFINEDSTATE1(X,T,NHAT,Q,THERMODYNAMICS_,    &
     &DIMENSIONLESS_,REFVALUES_)
              USE FLUIDDATA
              REAL(KIND=8), INTENT(IN) :: X(3)
              REAL(KIND=8), INTENT(IN) :: T
              REAL(KIND=8), INTENT(IN) :: NHAT(3)
              REAL(KIND=8), INTENT(INOUT) :: Q(5)
              TYPE (THERMODYNAMICS_T), INTENT(IN) :: THERMODYNAMICS_
              TYPE (DIMENSIONLESS_T), INTENT(IN) :: DIMENSIONLESS_
              TYPE (REFVALUES_T), INTENT(IN) :: REFVALUES_
            END SUBROUTINE USERDEFINEDSTATE1
          END INTERFACE 
        END MODULE USERDEFINEDSTATE1__genmod
