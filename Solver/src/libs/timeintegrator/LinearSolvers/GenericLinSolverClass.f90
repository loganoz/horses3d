!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!      GenericLinSolverClass.f90
!      Created: 2017-04-10 10:006:00 +0100 
!      By: Andrés Rueda
!
!      Class for defining common variables and type-bound procedures of linear solvers
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MODULE GenericLinSolverClass
   USE SMConstants
   USE DGSEMClass
   USE FTValueDictionaryClass
   use TimeIntegratorDefinitions
   IMPLICIT NONE
   
   PRIVATE
   PUBLIC GenericLinSolver_t, FTValueDictionary
   
   TYPE :: GenericLinSolver_t
      LOGICAL                                     :: converged = .FALSE.   ! The solution converged?
      INTEGER                                     :: DimPrb                ! Dimension of the problem
      INTEGER                                     :: niter = 0             ! Number of iterations to reach solution (for iterative solvers)
   CONTAINS
      !Subroutines:
      PROCEDURE :: construct
      PROCEDURE :: SetBValue
      PROCEDURE :: SetBValues
      PROCEDURE :: solve
      PROCEDURE :: GetXValue
      PROCEDURE :: destroy
      PROCEDURE :: SetOperatorDt
      PROCEDURE :: ReSetOperatorDt
      PROCEDURE :: AssemblyB
      !Functions:
      PROCEDURE :: Getxnorm    !Get solution norm
      PROCEDURE :: Getrnorm    !Get residual norm
      PROCEDURE :: ComputeANextStep
   END TYPE

CONTAINS

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE Construct(this,DimPrb,controlVariables,sem)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT), TARGET :: this
      INTEGER                  , INTENT(IN)            :: DimPrb
      TYPE(FTValueDictionary)  , INTENT(IN), OPTIONAL  :: controlVariables
      TYPE(DGSem), TARGET                  , OPTIONAL  :: sem
   END SUBROUTINE Construct
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE SetBValue(this, irow, value)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      INTEGER                  , INTENT(IN)  :: irow
      REAL(KIND=RP)            , INTENT(IN)  :: value
   END SUBROUTINE SetBValue

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE SetBValues(this, nvalues, irow, values)
      CLASS(GenericLinSolver_t)  , INTENT(INOUT)     :: this
      INTEGER                    , INTENT(IN)        :: nvalues
      INTEGER      , DIMENSION(:), INTENT(IN)        :: irow
      REAL(KIND=RP), DIMENSION(:), INTENT(IN)        :: values
   END SUBROUTINE SetBValues
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE solve(this,nEqn,nGradEqn,tol,maxiter,time,dt, ComputeTimeDerivative,computeA)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      integer,       intent(in)                :: nEqn
      integer,       intent(in)                :: nGradEqn
      REAL(KIND=RP), OPTIONAL                  :: tol
      INTEGER      , OPTIONAL                  :: maxiter
      REAL(KIND=RP), OPTIONAL                  :: time
      REAL(KIND=RP), OPTIONAL                  :: dt
      procedure(ComputeQDot_FCN)               :: ComputeTimeDerivative
      logical      , optional  , intent(inout) :: computeA
   END SUBROUTINE solve

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE GetXValue(this,irow,x_i)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      INTEGER                  , INTENT(IN)    :: irow
      REAL(KIND=RP)            , INTENT(OUT)   :: x_i
   END SUBROUTINE GetXValue

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE destroy(this)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
   END SUBROUTINE destroy

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE SetOperatorDt(this, dt)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)            , INTENT(IN)    :: dt
   END SUBROUTINE SetOperatorDt

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE ReSetOperatorDt(this, dt)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)            , INTENT(IN)    :: dt
   END SUBROUTINE ReSetOperatorDt
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE AssemblyB(this)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
   END SUBROUTINE AssemblyB
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   FUNCTION Getxnorm(this,TypeOfNorm) RESULT(xnorm)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      CHARACTER(len=*)                         :: TypeOfNorm
      REAL(KIND=RP)                            :: xnorm
   END FUNCTION Getxnorm
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   FUNCTION Getrnorm(this) RESULT(rnorm)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)                            :: rnorm
   END FUNCTION Getrnorm
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   FUNCTION ComputeANextStep(this) RESULT(ComputeA)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(IN) :: this
      LOGICAL                               :: ComputeA
   END FUNCTION ComputeANextStep
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
END MODULE GenericLinSolverClass
