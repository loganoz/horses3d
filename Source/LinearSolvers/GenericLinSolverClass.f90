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
   IMPLICIT NONE
   
   PRIVATE
   PUBLIC GenericLinSolver_t
   
   TYPE :: GenericLinSolver_t
      LOGICAL                                     :: converged = .FALSE.   ! The solution converged?
      INTEGER                                     :: DimPrb                ! Dimension of the problem
      INTEGER                                     :: niter = 0             ! Number of iterations to reach solution (for iterative solvers)
   CONTAINS
      !Subroutines:
      PROCEDURE :: construct
      PROCEDURE :: PreallocateA
      PROCEDURE :: ResetA
      PROCEDURE :: SetAColumn
      PROCEDURE :: AssemblyA
      PROCEDURE :: SetBValue
      PROCEDURE :: solve
      PROCEDURE :: GetCSRMatrix
      PROCEDURE :: GetXValue
      PROCEDURE :: destroy
      PROCEDURE :: SetOperatorDt
      PROCEDURE :: ReSetOperatorDt
      PROCEDURE :: AssemblyB
      !Functions:
      PROCEDURE :: Getxnorm    !Get solution norm
      PROCEDURE :: Getrnorm    !Get residual norm
   END TYPE

CONTAINS

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE Construct(this,DimPrb)
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      INTEGER                  , INTENT(IN)    :: DimPrb
   END SUBROUTINE Construct
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE PreallocateA(this,nnz)
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      INTEGER                                  :: nnz
   END SUBROUTINE PreallocateA
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE ResetA(this)
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
   END SUBROUTINE ResetA
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE SetAColumn(this,nvalues,irow,icol,values)
      CLASS(GenericLinSolver_t), INTENT(INOUT)  :: this
      INTEGER       , INTENT(IN)                :: nvalues
      INTEGER       , INTENT(IN), DIMENSION(:)  :: irow
      INTEGER       , INTENT(IN)                :: icol
      REAL(KIND=RP) , INTENT(IN), DIMENSION(:)  :: values
   END SUBROUTINE SetAColumn

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE AssemblyA(this)
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
   END SUBROUTINE AssemblyA

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE SetBValue(this, irow, value)
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      INTEGER                  , INTENT(IN)  :: irow
      REAL(KIND=RP)            , INTENT(IN)  :: value
   END SUBROUTINE SetBValue

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE solve(this,tol,maxiter)     
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP), OPTIONAL                  :: tol
      INTEGER      , OPTIONAL                  :: maxiter
   END SUBROUTINE solve

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE GetCSRMatrix(this,Acsr)         
      USE CSR_Matrices
      CLASS(GenericLinSolver_t), INTENT(IN)  :: this
      TYPE(csrMat)             , INTENT(OUT) :: Acsr 
   END SUBROUTINE GetCSRMatrix

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE GetXValue(this,irow,x_i)       
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      INTEGER                  , INTENT(IN)    :: irow
      REAL(KIND=RP)            , INTENT(OUT)   :: x_i
   END SUBROUTINE GetXValue

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE destroy(this)         
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
   END SUBROUTINE destroy

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE SetOperatorDt(this, dt)
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)            , INTENT(IN)    :: dt
   END SUBROUTINE SetOperatorDt

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE ReSetOperatorDt(this, dt)
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)            , INTENT(IN)    :: dt
   END SUBROUTINE ReSetOperatorDt
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   SUBROUTINE AssemblyB(this)         
       CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
   END SUBROUTINE AssemblyB
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   FUNCTION Getxnorm(this,TypeOfNorm) RESULT(xnorm)
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      CHARACTER(len=*)                         :: TypeOfNorm
      REAL(KIND=RP)                            :: xnorm
   END FUNCTION Getxnorm
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   FUNCTION Getrnorm(this) RESULT(rnorm)
      CLASS(GenericLinSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)                            :: rnorm
   END FUNCTION Getrnorm
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
END MODULE GenericLinSolverClass
