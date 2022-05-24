        !COMPILER-GENERATED INTERFACE MODULE: Tue May 24 12:21:52 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE USERDEFINEDNEUMANN1__genmod
          INTERFACE 
            SUBROUTINE USERDEFINEDNEUMANN1(X,T,NHAT,Q,U_X,U_Y,U_Z,FLUX, &
     &THERMODYNAMICS_,DIMENSIONLESS_,REFVALUES_)
              USE FLUIDDATA
              REAL(KIND=8), INTENT(IN) :: X(3)
              REAL(KIND=8), INTENT(IN) :: T
              REAL(KIND=8), INTENT(IN) :: NHAT(3)
              REAL(KIND=8), INTENT(IN) :: Q(5)
              REAL(KIND=8), INTENT(IN) :: U_X(5)
              REAL(KIND=8), INTENT(IN) :: U_Y(5)
              REAL(KIND=8), INTENT(IN) :: U_Z(5)
              REAL(KIND=8), INTENT(INOUT) :: FLUX(5)
              TYPE (THERMODYNAMICS_T), INTENT(IN) :: THERMODYNAMICS_
              TYPE (DIMENSIONLESS_T), INTENT(IN) :: DIMENSIONLESS_
              TYPE (REFVALUES_T), INTENT(IN) :: REFVALUES_
            END SUBROUTINE USERDEFINEDNEUMANN1
          END INTERFACE 
        END MODULE USERDEFINEDNEUMANN1__genmod
