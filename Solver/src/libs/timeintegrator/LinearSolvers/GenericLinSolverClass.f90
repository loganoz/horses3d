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
   PUBLIC GenericLinSolver_t, FTValueDictionary, MatrixShift_FCN, Default_MatrixShift, MatrixShift
   
   TYPE :: GenericLinSolver_t
      LOGICAL                                     :: converged = .FALSE.   ! The solution converged?
      INTEGER                                     :: DimPrb                ! Dimension of the problem
      INTEGER                                     :: niter = 0             ! Number of iterations to reach solution (for iterative solvers)
   CONTAINS
      !Subroutines:
      PROCEDURE :: construct
      PROCEDURE :: SetRHSValue
      PROCEDURE :: SetRHSValues
      PROCEDURE :: SetRHS
      PROCEDURE :: solve
      PROCEDURE :: GetXValue
      PROCEDURE :: GetX
      PROCEDURE :: destroy
      PROCEDURE :: SetOperatorDt
      PROCEDURE :: ReSetOperatorDt
      PROCEDURE :: AssemblyRHS
      !Functions:
      PROCEDURE :: Getxnorm    !Get solution norm
      PROCEDURE :: Getrnorm    !Get residual norm
      PROCEDURE :: ComputeANextStep
   END TYPE
   
   abstract interface
      function MatrixShift_FCN(dt) result(Ashift)
         use SMConstants
         implicit none
         !------------------------------------------------------
         real(kind=RP), intent(in) :: dt
         real(kind=RP)             :: Ashift
         !------------------------------------------------------
      end function MatrixShift_FCN
   end interface
   
   procedure(MatrixShift_FCN), pointer :: MatrixShift

CONTAINS
   
   function Default_MatrixShift(dt) result(Ashift)
      use SMConstants
      implicit none
      !------------------------------------------------------
      real(kind=RP), intent(in) :: dt
      real(kind=RP)             :: Ashift
      !------------------------------------------------------
      print*, 'using default'
      ! Do nothing
      Ashift = 0._RP
   end function Default_MatrixShift
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE Construct(this,DimPrb,controlVariables,sem,MatrixShiftFunc)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT), TARGET :: this
      INTEGER                  , INTENT(IN)            :: DimPrb
      TYPE(FTValueDictionary)  , INTENT(IN), OPTIONAL  :: controlVariables
      TYPE(DGSem), TARGET                  , OPTIONAL  :: sem
      procedure(MatrixShift_FCN)                       :: MatrixShiftFunc
      
      ERROR stop ':: Linear solver does not have a constructor yet'
   END SUBROUTINE Construct
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE SetRHS(this, RHS)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)            , INTENT(IN)    :: RHS(:)
      
      ERROR stop ':: SetRHS not implemented for desired linear solver'
   END SUBROUTINE SetRHS

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE SetRHSValue(this, irow, value)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      INTEGER                  , INTENT(IN)  :: irow
      REAL(KIND=RP)            , INTENT(IN)  :: value
      
      ERROR stop ':: SetRHSValue not implemented for desired linear solver'
   END SUBROUTINE SetRHSValue

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE SetRHSValues(this, nvalues, irow, values)
      CLASS(GenericLinSolver_t)  , INTENT(INOUT)     :: this
      INTEGER                    , INTENT(IN)        :: nvalues
      INTEGER      , DIMENSION(:), INTENT(IN)        :: irow
      REAL(KIND=RP), DIMENSION(:), INTENT(IN)        :: values
      
      ERROR stop ':: SetRHSValues not implemented for desired linear solver'
   END SUBROUTINE SetRHSValues
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE solve(this,nEqn, nGradEqn, ComputeTimeDerivative,tol,maxiter,time,dt,computeA)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      integer,       intent(in)                :: nEqn
      integer,       intent(in)                :: nGradEqn
      procedure(ComputeQDot_FCN)               :: ComputeTimeDerivative
      REAL(KIND=RP), OPTIONAL                  :: tol
      INTEGER      , OPTIONAL                  :: maxiter
      REAL(KIND=RP), OPTIONAL                  :: time
      REAL(KIND=RP), OPTIONAL                  :: dt
      logical      , optional  , intent(inout) :: computeA
      
      ERROR stop ':: solve not implemented for desired linear solver!!!'
   END SUBROUTINE solve
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE GetXValue(this,irow,x_i)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      INTEGER                  , INTENT(IN)    :: irow
      REAL(KIND=RP)            , INTENT(OUT)   :: x_i
      
      ERROR stop ':: GetXValue not implemented for desired linear solver'
   END SUBROUTINE GetXValue

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   function GetX(this) result(x)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)                            :: x(this % DimPrb)
      
      ERROR stop ':: GetX not implemented for desired linear solver'
   end function GetX

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE destroy(this)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      
      write(STD_OUT,*) 'WARNING :: destroy not implemented for desired linear solver'
   END SUBROUTINE destroy

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE SetOperatorDt(this, dt)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)            , INTENT(IN)    :: dt
      
      write(STD_OUT,*) 'WARNING :: SetOperatorDt not implemented for desired linear solver'
   END SUBROUTINE SetOperatorDt

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE ReSetOperatorDt(this, dt)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)            , INTENT(IN)    :: dt
      
      write(STD_OUT,*) 'WARNING :: ReSetOperatorDt not implemented for desired linear solver'
   END SUBROUTINE ReSetOperatorDt
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE AssemblyRHS(this)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
   END SUBROUTINE AssemblyRHS
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   FUNCTION Getxnorm(this,TypeOfNorm) RESULT(xnorm)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      CHARACTER(len=*)                         :: TypeOfNorm
      REAL(KIND=RP)                            :: xnorm
      
      ERROR stop ':: Getxnorm not implemented for desired linear solver'
   END FUNCTION Getxnorm
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   FUNCTION Getrnorm(this) RESULT(rnorm)
      IMPLICIT NONE
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)                            :: rnorm
      
      ERROR stop ':: Getrnorm not implemented for desired linear solver'
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
