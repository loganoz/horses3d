        !COMPILER-GENERATED INTERFACE MODULE: Tue May 17 10:06:14 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE USERDEFINEDSOURCETERMNS__genmod
          INTERFACE 
            SUBROUTINE USERDEFINEDSOURCETERMNS(X,Q,TIME,S,              &
     &THERMODYNAMICS_,DIMENSIONLESS_,REFVALUES_)
              USE FLUIDDATA
              REAL(KIND=8), INTENT(IN) :: X(3)
              REAL(KIND=8), INTENT(IN) :: Q(5)
              REAL(KIND=8), INTENT(IN) :: TIME
              REAL(KIND=8), INTENT(INOUT) :: S(5)
              TYPE (THERMODYNAMICS_T), INTENT(IN) :: THERMODYNAMICS_
              TYPE (DIMENSIONLESS_T), INTENT(IN) :: DIMENSIONLESS_
              TYPE (REFVALUES_T), INTENT(IN) :: REFVALUES_
            END SUBROUTINE USERDEFINEDSOURCETERMNS
          END INTERFACE 
        END MODULE USERDEFINEDSOURCETERMNS__genmod
