        !COMPILER-GENERATED INTERFACE MODULE: Tue May 17 10:06:12 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE USERDEFINEDGRADVARS1__genmod
          INTERFACE 
            SUBROUTINE USERDEFINEDGRADVARS1(X,T,NHAT,Q,U,GETGRADIENTS,  &
     &THERMODYNAMICS_,DIMENSIONLESS_,REFVALUES_)
              USE FLUIDDATA
              REAL(KIND=8), INTENT(IN) :: X(3)
              REAL(KIND=8), INTENT(IN) :: T
              REAL(KIND=8), INTENT(IN) :: NHAT(3)
              REAL(KIND=8), INTENT(IN) :: Q(5)
              REAL(KIND=8), INTENT(INOUT) :: U(5)
              INTERFACE 
                SUBROUTINE GETGRADIENTS(NEQN,NGRADEQN,Q,U,RHO_)
                  USE FLUIDDATA
                  INTEGER(KIND=4), INTENT(IN) :: NEQN
                  INTEGER(KIND=4), INTENT(IN) :: NGRADEQN
                  REAL(KIND=8), INTENT(IN) :: Q(NEQN)
                  REAL(KIND=8), INTENT(OUT) :: U(NGRADEQN)
                  REAL(KIND=8) ,OPTIONAL, INTENT(IN) :: RHO_
                END SUBROUTINE GETGRADIENTS
              END INTERFACE 
              TYPE (THERMODYNAMICS_T), INTENT(IN) :: THERMODYNAMICS_
              TYPE (DIMENSIONLESS_T), INTENT(IN) :: DIMENSIONLESS_
              TYPE (REFVALUES_T), INTENT(IN) :: REFVALUES_
            END SUBROUTINE USERDEFINEDGRADVARS1
          END INTERFACE 
        END MODULE USERDEFINEDGRADVARS1__genmod
