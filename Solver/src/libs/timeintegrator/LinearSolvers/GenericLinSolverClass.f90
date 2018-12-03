!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!      GenericLinSolverClass.f90
!      Created: 2017-04-10 10:006:00 +0100 
!      By: Andrés Rueda
!
!      Class for defining common variables and type-bound procedures of linear solvers
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
module GenericLinSolverClass
   use SMConstants
   use DGSEMClass
   use FTValueDictionaryClass
   use TimeIntegratorDefinitions
   use MatrixClass         , only: Matrix_t
   use AnalyticalJacobian  , only: AnalyticalJacobian_Compute
   use NumericalJacobian   , only: NumericalJacobian_Compute
   implicit none
   
   private
   public GenericLinSolver_t
   public MatrixShift_FCN
   public Default_MatrixShift, MatrixShift
   public NUMERICAL_JACOBIAN, ANALYTICAL_JACOBIAN
   
   public FTValueDictionary
   
   integer, parameter :: NUMERICAL_JACOBIAN  = 1
   integer, parameter :: ANALYTICAL_JACOBIAN = 2
   
   type :: GenericLinSolver_t
      logical              :: converged = .FALSE.   ! The solution converged?
      integer              :: DimPrb                ! Dimension of the problem
      integer              :: niter = 0             ! Number of iterations to reach solution (for iterative solvers)
      integer              :: JacobianComputation = NUMERICAL_JACOBIAN
      type(DGSem), pointer :: p_sem   
   contains
      !Subroutines:
      procedure :: construct
      procedure :: SetRHSValue
      procedure :: SetRHSValues
      procedure :: SetRHS
      procedure :: solve
      procedure :: GetXValue
      procedure :: GetX
      procedure :: destroy
      procedure :: SetOperatorDt
      procedure :: ReSetOperatorDt
      procedure :: AssemblyRHS
      procedure :: ComputeJacobian
      !Functions:
      procedure :: Getxnorm    !Get solution norm
      procedure :: Getrnorm    !Get residual norm
      procedure :: ComputeANextStep
   end type
   
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
   
   procedure(MatrixShift_FCN), pointer :: MatrixShift   ! TODO?: move to GenericLinSolver_t to allow different MatrixShifts for different solvers?

contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
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
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Construct(this,DimPrb,controlVariables,sem,MatrixShiftFunc)
      implicit none
      !-arguments-----------------------------------------------------------
      class(GenericLinSolver_t), intent(inout), target :: this
      integer                  , intent(in)            :: DimPrb
      type(FTValueDictionary)  , intent(in), optional  :: controlVariables
      type(DGSem), target                  , optional  :: sem
      procedure(MatrixShift_FCN)                       :: MatrixShiftFunc
      !---------------------------------------------------------------------
      ERROR stop ':: Linear solver does not have a constructor yet'
   end subroutine Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetRHS(this, RHS)
      implicit none
      class(GenericLinSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: RHS(:)
      
      ERROR stop ':: SetRHS not implemented for desired linear solver'
   end subroutine SetRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetRHSValue(this, irow, value)
      implicit none
      class(GenericLinSolver_t), intent(inout) :: this
      integer                  , intent(in)  :: irow
      real(kind=RP)            , intent(in)  :: value
      
      ERROR stop ':: SetRHSValue not implemented for desired linear solver'
   end subroutine SetRHSValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetRHSValues(this, nvalues, irow, values)
      class(GenericLinSolver_t)  , intent(inout)     :: this
      integer                    , intent(in)        :: nvalues
      integer      , DIMENSION(:), intent(in)        :: irow
      real(kind=RP), DIMENSION(:), intent(in)        :: values
      
      ERROR stop ':: SetRHSValues not implemented for desired linear solver'
   end subroutine SetRHSValues
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine solve(this,nEqn, nGradEqn, ComputeTimeDerivative,tol,maxiter,time,dt,computeA)
      implicit none
      class(GenericLinSolver_t), intent(inout) :: this
      integer,       intent(in)                :: nEqn
      integer,       intent(in)                :: nGradEqn
      procedure(ComputeTimeDerivative_f)               :: ComputeTimeDerivative
      real(kind=RP), optional                  :: tol
      integer      , optional                  :: maxiter
      real(kind=RP), optional                  :: time
      real(kind=RP), optional                  :: dt
      logical      , optional  , intent(inout) :: computeA
      
      ERROR stop ':: solve not implemented for desired linear solver!!!'
   end subroutine solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine GetXValue(this,irow,x_i)
      implicit none
      class(GenericLinSolver_t), intent(inout) :: this
      integer                  , intent(in)    :: irow
      real(kind=RP)            , intent(OUT)   :: x_i
      
      ERROR stop ':: GetXValue not implemented for desired linear solver'
   end subroutine GetXValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function GetX(this) result(x)
      implicit none
      class(GenericLinSolver_t), intent(inout) :: this
      real(kind=RP)                            :: x(this % DimPrb)
      
      ERROR stop ':: GetX not implemented for desired linear solver'
   end function GetX
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine destroy(this)
      implicit none
      class(GenericLinSolver_t), intent(inout) :: this
      
      write(STD_OUT,*) 'WARNING :: destroy not implemented for desired linear solver'
   end subroutine destroy
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetOperatorDt(this, dt)
      implicit none
      class(GenericLinSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: dt
      
      write(STD_OUT,*) 'WARNING :: SetOperatorDt not implemented for desired linear solver'
   end subroutine SetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ReSetOperatorDt(this, dt)
      implicit none
      class(GenericLinSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: dt
      
      write(STD_OUT,*) 'WARNING :: ReSetOperatorDt not implemented for desired linear solver'
   end subroutine ReSetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine AssemblyRHS(this)
      implicit none
      class(GenericLinSolver_t), intent(inout) :: this
   end subroutine AssemblyRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function Getxnorm(this,TypeOfNorm) RESULT(xnorm)
      implicit none
      class(GenericLinSolver_t), intent(inout) :: this
      character(len=*)                         :: TypeOfNorm
      real(kind=RP)                            :: xnorm
      
      ERROR stop ':: Getxnorm not implemented for desired linear solver'
   end function Getxnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function Getrnorm(this) RESULT(rnorm)
      implicit none
      class(GenericLinSolver_t), intent(inout) :: this
      real(kind=RP)                            :: rnorm
      
      ERROR stop ':: Getrnorm not implemented for desired linear solver'
   end function Getrnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function ComputeANextStep(this) RESULT(ComputeA)
      implicit none
      class(GenericLinSolver_t), intent(in) :: this
      logical                               :: ComputeA
   end function ComputeANextStep
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ComputeJacobian(this,Matrix,dt,time,nEqn,nGradEqn,ComputeTimeDerivative)
      implicit none
      !-----------------------------------------------------------
      class(GenericLinSolver_t), intent(inout) :: this
      class(Matrix_t)                          :: Matrix
      real(kind=RP), intent(in)                :: dt
      real(kind=RP), intent(in)                :: time
      integer,       intent(in)                :: nEqn
      integer,       intent(in)                :: nGradEqn
      procedure(ComputeTimeDerivative_f)       :: ComputeTimeDerivative
      !-----------------------------------------------------------
      
      select case (this % JacobianComputation)
         case(NUMERICAL_JACOBIAN)
            call NumericalJacobian_Compute(this % p_sem, nEqn, nGradEqn, time, Matrix, ComputeTimeDerivative, .TRUE. )
         case(ANALYTICAL_JACOBIAN)
            call AnalyticalJacobian_Compute(this % p_sem, nEqn, time, Matrix)
      end select
         
      call Matrix % shift( MatrixShift(dt) )
      
   end subroutine ComputeJacobian
end module GenericLinSolverClass
