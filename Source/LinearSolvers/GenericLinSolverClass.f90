!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!      GenericLinSolverClass.f90
!      Created: 2017-04-10 10:006:00 +0100 
!      By: Andrés Rueda
!
!      Class for defining type-bound procedures of linear solvers
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MODULE GenericLinSolverClass
   USE SMConstants
   IMPLICIT NONE
   
   TYPE, ABSTRACT :: GenericLinSolver_t
      LOGICAL                                     :: converged = .FALSE.   ! The solution converged?
      INTEGER                                     :: DimPrb                ! Dimension of the problem
      INTEGER                                     :: niter                 ! Number of iterations to reach solution (for iterative solvers)
   CONTAINS
      !Subroutines:
      PROCEDURE(IConstruct)      , pass, DEFERRED :: construct
      PROCEDURE(IPreallocateA)   , pass, DEFERRED :: PreallocateA
      PROCEDURE(IResetA)         , pass, DEFERRED :: ResetA
      PROCEDURE(ISetAColumn)     , pass, DEFERRED :: SetAColumn
      PROCEDURE(IAssemblyA)      , pass, DEFERRED :: AssemblyA
      PROCEDURE(ISetBValue)      , pass, DEFERRED :: SetBValue
      PROCEDURE(Isolve)          , pass, DEFERRED :: solve
      PROCEDURE(IGetCSRMatrix)   , pass, DEFERRED :: GetCSRMatrix
      PROCEDURE(IGetXValue)      , pass, DEFERRED :: GetXValue
      PROCEDURE(Idestruct)       , pass, DEFERRED :: destroy
      PROCEDURE(ISetOperatorDt)  , pass, DEFERRED :: SetOperatorDt
      PROCEDURE(IReSetOperatorDt), pass, DEFERRED :: ReSetOperatorDt
      PROCEDURE(IAssemblyB)      , pass, DEFERRED :: AssemblyB
      !Functions:
      PROCEDURE(IGetxnorm)       , pass, DEFERRED :: Getxnorm    !Get solution norm
      PROCEDURE(IGetrnorm)       , pass, DEFERRED :: Getrnorm    !Get residual norm
   END TYPE
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Interface for Construct procedure
   abstract interface
      SUBROUTINE IConstruct(this,DimPrb)         
         import :: GenericLinSolver_t
         class(GenericLinSolver_t), INTENT(INOUT) :: this
         INTEGER                  , INTENT(IN)    :: DimPrb
      end SUBROUTINE IConstruct
   end interface
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Interface for PreallocateA procedure
   abstract interface
      SUBROUTINE IPreallocateA(this,nnz)         
         import :: GenericLinSolver_t
         class(GenericLinSolver_t), INTENT(INOUT) :: this
         INTEGER                                  :: nnz
      end SUBROUTINE IPreallocateA
   end interface
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Interface for ResetA procedure
   abstract interface
      SUBROUTINE IResetA(this)         
         import :: GenericLinSolver_t
         class(GenericLinSolver_t), INTENT(INOUT) :: this
      end SUBROUTINE IResetA
   end interface
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Interface for SetAColumn procedure
   abstract interface
      SUBROUTINE ISetAColumn(this,nvalues,irow,icol,values)         
         USE SMConstants
         import :: GenericLinSolver_t
         
         class(GenericLinSolver_t), INTENT(INOUT)  :: this
         INTEGER       , INTENT(IN)                :: nvalues
         INTEGER       , INTENT(IN), DIMENSION(:)  :: irow
         INTEGER       , INTENT(IN)                :: icol
         REAL(KIND=RP) , INTENT(IN), DIMENSION(:)  :: values
      end SUBROUTINE ISetAColumn
   end interface
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Interface for AssemblyA procedure
   abstract interface
      SUBROUTINE IAssemblyA(this)         
         import :: GenericLinSolver_t
         class(GenericLinSolver_t), INTENT(INOUT) :: this
      end SUBROUTINE IAssemblyA
   end interface
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Interface for SetBValue procedure
   abstract interface
      SUBROUTINE ISetBValue(this, irow, value)
         USE SMConstants
         import :: GenericLinSolver_t
         
         class(GenericLinSolver_t), INTENT(INOUT) :: this
         INTEGER                  , INTENT(IN)  :: irow
         REAL(KIND=RP)            , INTENT(IN)  :: value
      end SUBROUTINE ISetBValue
   end interface
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Interface for solve procedure
   abstract interface
      SUBROUTINE Isolve(this,tol,maxiter)         
         USE SMConstants
         import :: GenericLinSolver_t
         
         class(GenericLinSolver_t), INTENT(INOUT) :: this
         REAL(KIND=RP), OPTIONAL                  :: tol
         INTEGER      , OPTIONAL                  :: maxiter
      end SUBROUTINE Isolve
   end interface
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Interface for GetCSRMatrix procedure
   abstract interface
      SUBROUTINE IGetCSRMatrix(this,Acsr)         
         USE CSR_Matrices
         import :: GenericLinSolver_t
         
         class(GenericLinSolver_t), INTENT(IN)  :: this
         TYPE(csrMat)             , INTENT(OUT) :: Acsr 
      end SUBROUTINE IGetCSRMatrix
   end interface
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Interface for GetXValue procedure
   abstract interface
      SUBROUTINE IGetXValue(this,irow,x_i)       
         USE SMConstants
         import :: GenericLinSolver_t
         
         class(GenericLinSolver_t), INTENT(INOUT) :: this
         INTEGER                  , INTENT(IN)    :: irow
         REAL(KIND=RP)            , INTENT(OUT)   :: x_i
      end SUBROUTINE IGetXValue
   end interface
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Interface for destruct procedure
   abstract interface
      SUBROUTINE Idestruct(this)         
         import :: GenericLinSolver_t
         class(GenericLinSolver_t), INTENT(INOUT) :: this
      end SUBROUTINE Idestruct
   end interface
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Interface for SetOperatorDt procedure
   abstract interface
      SUBROUTINE ISetOperatorDt(this, dt)
         USE SMConstants
         import :: GenericLinSolver_t
         
         class(GenericLinSolver_t), INTENT(INOUT) :: this
         REAL(KIND=RP)            , INTENT(IN)    :: dt
      end SUBROUTINE ISetOperatorDt
   end interface
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Interface for ReSetOperatorDt procedure
   abstract interface
      SUBROUTINE IReSetOperatorDt(this, dt)
         USE SMConstants
         import :: GenericLinSolver_t
         
         class(GenericLinSolver_t), INTENT(INOUT) :: this
         REAL(KIND=RP)            , INTENT(IN)    :: dt
      end SUBROUTINE IReSetOperatorDt
   end interface
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Interface for AssemblyB procedure
   abstract interface
      SUBROUTINE IAssemblyB(this)         
         import :: GenericLinSolver_t
         class(GenericLinSolver_t), INTENT(INOUT) :: this
      end SUBROUTINE IAssemblyB
   end interface
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Interface for Getxnorm procedure
   abstract interface
      FUNCTION IGetxnorm(this,TypeOfNorm) RESULT(xnorm)
         USE SMConstants
         import :: GenericLinSolver_t
         
         class(GenericLinSolver_t), INTENT(INOUT) :: this
         CHARACTER(len=*)                         :: TypeOfNorm
         REAL(KIND=RP)                            :: xnorm
      END FUNCTION IGetxnorm
   end interface
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Interface for Getrnorm procedure
   abstract interface
      FUNCTION IGetrnorm(this) RESULT(rnorm)
         USE SMConstants
         import :: GenericLinSolver_t
         
         class(GenericLinSolver_t), INTENT(INOUT) :: this
         REAL(KIND=RP)                            :: rnorm
      END FUNCTION IGetrnorm
   end interface
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
END MODULE GenericLinSolverClass
