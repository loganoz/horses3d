!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!      Class for defining common variables and type-bound procedures of linear solvers
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
module GenericLinSolverClass
   use SMConstants
   use DGSEMClass
   use FTValueDictionaryClass
   use TimeIntegratorDefinitions
   use MatrixClass            , only: Matrix_t
   use AnalyticalJacobian     , only: AnJacobian_t
   use NumericalJacobian      , only: NumJacobian_t
   use JacobianComputerClass  , only: JacobianComputer_t, GetJacobianFlag
   use MPI_Process_Info       , only: MPI_Process
   implicit none
   
   private
   public GenericLinSolver_t
   public MatrixShift_FCN
   public Default_MatrixShift, MatrixShift
   public NUMERICAL_JACOBIAN, ANALYTICAL_JACOBIAN
   
   public FTValueDictionary
   
   
   type :: GenericLinSolver_t
      class(JacobianComputer_t), allocatable   :: Jacobian
      logical                          :: converged = .FALSE.   ! The solution converged?
      logical                          :: withMPI = .FALSE.
      integer                          :: DimPrb                ! Dimension of the (local) problem
      integer                          :: globalDimPrb          ! Dimension of the (global) problem
      integer                          :: niter = 0             ! Number of iterations to reach solution (for iterative solvers)
      integer                          :: JacobianComputation = NUMERICAL_JACOBIAN
      type(DGSem), pointer             :: p_sem => null()
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
      procedure :: SetJacobian
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
   
   procedure(MatrixShift_FCN), pointer :: MatrixShift =>  Default_MatrixShift  ! TODO?: move to GenericLinSolver_t to allow different MatrixShifts for different solvers?

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
      
      ! Do nothing
      Ashift = 0._RP
   end function Default_MatrixShift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Construct(this, DimPrb, globalDimPrb, nEqn, controlVariables, sem, MatrixShiftFunc)
      implicit none
      !-arguments-----------------------------------------------------------
      class(GenericLinSolver_t), intent(inout), target :: this
      integer                  , intent(in)            :: DimPrb
      integer                  , intent(in)            :: globalDimPrb        
      integer                  , intent(in)            :: nEqn
      type(FTValueDictionary)  , intent(in), optional  :: controlVariables
      type(DGSem), target                  , optional  :: sem
      procedure(MatrixShift_FCN)                       :: MatrixShiftFunc     ! TODO: Make this optional
      !---------------------------------------------------------------------
      
      if (globalDimPrb < DimPrb) then        ! This never makes sense
         error stop 'Inconsistent problem sizes: globalDimPrb < DimPrb'
      elseif (globalDimPrb > DimPrb) then    ! This only makes sense if MPI is active
         if (.not. MPI_Process % doMPIAction) then    ! MPI is not enabled: ERROR
            error stop "Trying to solve linSystem with MPI, but there's no MPI"
         end if
         this % withMPI = .TRUE.
      end if
      
      this % JacobianComputation = GetJacobianFlag()
      
      select case (this % JacobianComputation)
         case (NOTDEF_JACOBIAN )    ; allocate(this % Jacobian)
         case (NUMERICAL_JACOBIAN ) ; allocate(NumJacobian_t :: this % Jacobian)
         case (ANALYTICAL_JACOBIAN) ; allocate(AnJacobian_t  :: this % Jacobian)
         case default
            error stop 'Invalid jacobian type'
      end select
      
!
!     ***************************
!     Construct Jacobian computer
!     ***************************
!
      if ( present(sem) ) then
         call this % Jacobian % construct(sem % mesh, nEqn, controlVariables)
      end if
   end subroutine Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetRHS(this, RHS)
      implicit none
      class(GenericLinSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: RHS(this % DimPrb)
      
      error stop ':: SetRHS not implemented for desired linear solver'
   end subroutine SetRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetJacobian(this,Matrix)
      implicit none
      !-arguments-----------------------------------------------------------
      class(GenericLinSolver_t), intent(inout)  :: this
      class(Matrix_t)          , intent(in)     :: Matrix
      !---------------------------------------------------------------------
      
      error stop ':: SetJacobian not implemented for desired linear solver'
      
   end subroutine SetJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetRHSValue(this, irow, value)
      implicit none
      class(GenericLinSolver_t), intent(inout) :: this
      integer                  , intent(in)  :: irow
      real(kind=RP)            , intent(in)  :: value
      
      error stop ':: SetRHSValue not implemented for desired linear solver'
   end subroutine SetRHSValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetRHSValues(this, nvalues, irow, values)
      class(GenericLinSolver_t)  , intent(inout)     :: this
      integer                    , intent(in)        :: nvalues
      integer      , DIMENSION(:), intent(in)        :: irow
      real(kind=RP), DIMENSION(:), intent(in)        :: values
      
      error stop ':: SetRHSValues not implemented for desired linear solver'
   end subroutine SetRHSValues
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine solve(this,nEqn, nGradEqn, ComputeTimeDerivative,tol,maxiter,time,dt,computeA)
      implicit none
      class(GenericLinSolver_t), target, intent(inout) :: this
      integer,       intent(in)                :: nEqn
      integer,       intent(in)                :: nGradEqn
      procedure(ComputeTimeDerivative_f)       :: ComputeTimeDerivative
      real(kind=RP), optional                  :: tol
      integer      , optional                  :: maxiter
      real(kind=RP), optional                  :: time
      real(kind=RP), optional                  :: dt
      logical      , optional  , intent(inout) :: computeA
      
      error stop ':: solve not implemented for desired linear solver!!!'
   end subroutine solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine GetXValue(this,irow,x_i)
      implicit none
      class(GenericLinSolver_t), intent(inout) :: this
      integer                  , intent(in)    :: irow
      real(kind=RP)            , intent(OUT)   :: x_i
      
      error stop ':: GetXValue not implemented for desired linear solver'
   end subroutine GetXValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function GetX(this) result(x)
      implicit none
      class(GenericLinSolver_t), intent(inout) :: this
      real(kind=RP)                            :: x(this % DimPrb)
      
      error stop ':: GetX not implemented for desired linear solver'
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
      
      error stop ':: Getxnorm not implemented for desired linear solver'
   end function Getxnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function Getrnorm(this) RESULT(rnorm)
      implicit none
      class(GenericLinSolver_t), intent(inout) :: this
      real(kind=RP)                            :: rnorm
      
      error stop ':: Getrnorm not implemented for desired linear solver'
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
end module GenericLinSolverClass