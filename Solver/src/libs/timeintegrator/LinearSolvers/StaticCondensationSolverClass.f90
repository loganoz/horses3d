!
!//////////////////////////////////////////////////////
!
!   @File:    StaticCondensationSolverClass.f90
!   @Author:  Andrés Rueda (am.rueda@upm.es)
!   @Created: Tue Dec  4 16:26:02 2018
!   @Last revision date: Tue Dec  4 21:53:50 2018
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 9b3844379fde2350e64816efcdf3bf724c8b3041
!
!//////////////////////////////////////////////////////
!
!  StaticCondensationSolverClass: (UNDER DEVELOPMENT!!)
!     Routines for solving Gauss-Lobatto DGSEM representations using static-condensation/substructuring
!
!     -> Only valid for hp-conforming representations (a middle ground can be found -partial condensation)
!
!  TO CHANGE:
!     -> Uniform polynomial order is supposed in the moment
!     -> nEqn hard-coded to 5
!     -> Isolating DOFs on physical boundaries (can be avoided)
!
module StaticCondensationSolverClass
   use SMConstants
   use DGSEMClass                , only: DGSem, ComputeTimeDerivative_f
   use GenericLinSolverClass
   use StaticCondensedMatrixClass, only: StaticCondensedMatrix_t
   implicit none
   
   private
   public StaticCondSolver_t
   
   type, extends(GenericLinSolver_t) :: StaticCondSolver_t
      type(StaticCondensedMatrix_t) :: A
      real(kind=RP)                 :: Ashift = 0._RP ! Current shift of the Jacobian matrix
      real(kind=RP), allocatable    :: x(:)           ! Solution vector
      real(kind=RP), allocatable    :: b(:)           ! Right hand side
   contains
      procedure :: construct        => SCS_construct
      procedure :: destroy          => SCS_destruct
      procedure :: SetOperatorDt    => SCS_SetOperatorDt
      procedure :: ReSetOperatorDt  => SCS_ReSetOperatorDt
      procedure :: solve            => SCS_solve
      procedure :: SetRHS           => SCS_SetRHS
   end type StaticCondSolver_t
   
   integer, parameter :: nEqn = 5 ! hard-coded
   
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------
!  Constructor
!  -----------
   subroutine SCS_construct(this,DimPrb,controlVariables,sem,MatrixShiftFunc)
      implicit none
      !-arguments-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout), target :: this
      integer                  , intent(in)            :: DimPrb
      type(FTValueDictionary)  , intent(in), optional  :: controlVariables
      type(DGSem), target                  , optional  :: sem
      procedure(MatrixShift_FCN)                       :: MatrixShiftFunc
      !-local-variables-----------------------------------------------------
      integer :: size_i
      integer :: N ! Polynomial order
      integer :: nelem
      !---------------------------------------------------------------------
      
!
!     ----------------------
!     Check needed arguments
!     ----------------------
!
      if (.not. present(controlVariables)) ERROR stop 'StaticCondSolver_t needs controlVariables'
      if (.not. present(sem)) ERROR stop 'StaticCondSolver_t needs DGSem'
      
      ! TODO: Add conformity check!
      ! TODO: Add GL check
      
      this % DimPrb = DimPrb
      this % p_sem => sem
      
      if ( controlVariables % containsKey("jacobian flag") ) then
         this % JacobianComputation = controlVariables % integerValueForKey("jacobian flag")
      end if
      
      MatrixShift => MatrixShiftFunc
      
!
!     -----------
!     Allocations
!     -----------
!
      allocate(this % x(DimPrb))
      allocate(this % b(DimPrb))
      
      
      N = sem % mesh % elements(1) % Nxyz(1) ! hard-coded: Uniform order as elem 1-x
      nelem  = sem % mesh % no_of_elements
      size_i = nelem * (N-1)**3 * nEqn
      
      call this % A % construct  (num_of_Rows = DimPrb, &
                                  num_of_Blocks = nelem, &
                                  num_of_Rows_reduced = size_i )
!
!     -------------------------
!     Construct the permutation
!     -------------------------
!
      call this % A % constructPermutationArrays(sem % mesh % Nx, sem % mesh % Ny, sem % mesh % Nz, nEqn)
      
   end subroutine SCS_construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------
!  Destructor
!  ----------
   subroutine SCS_destruct(this)
      implicit none
      !-arguments-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout) :: this
      !---------------------------------------------------------------------
      
      call this % A % destruct
      
      deallocate (this % x)
      deallocate (this % b)
   end subroutine SCS_destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SCS_SetOperatorDt(this,dt)       
      implicit none
      !-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: dt
      !-----------------------------------------------------------
      
      this % Ashift = MatrixShift(dt)
      call this % A % Shift(this % Ashift)
      
    end subroutine SCS_SetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SCS_ReSetOperatorDt(this,dt)       
      implicit none
      !-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: dt
      !-----------------------------------------------------------
      real(kind=RP)                            :: shift
      !-----------------------------------------------------------
      
      shift = MatrixShift(dt)
      call this % A % Shift(-this % Ashift)
      call this % A % Shift(shift)
      this % Ashift = shift
      
   end subroutine SCS_ReSetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SCS_SetRHS(this, RHS)
      implicit none
      !-arguments-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: RHS(this % DimPrb)
      !---------------------------------------------------------------------
      
      this % b = RHS
      
   end subroutine SCS_SetRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SCS_solve(this,nEqn, nGradEqn, ComputeTimeDerivative,tol,maxiter,time,dt,computeA)
      use CSRMatrixClass ! debug
      implicit none
      !-arguments-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout) :: this
      integer,       intent(in)                :: nEqn
      integer,       intent(in)                :: nGradEqn
      procedure(ComputeTimeDerivative_f)       :: ComputeTimeDerivative
      real(kind=RP), optional                  :: tol
      integer      , optional                  :: maxiter
      real(kind=RP), optional                  :: time
      real(kind=RP), optional                  :: dt
      logical      , optional  , intent(inout) :: computeA
      !-local-variables-----------------------------------------------------
      type(csrMat_t) :: A ! debug
      !---------------------------------------------------------------------
      
!     Compute Jacobian matrix if needed
!     ---------------------------------      
      if ( present(ComputeA)) then
         if (ComputeA) then
            call this % ComputeJacobian(this % A,dt,time,nEqn,nGradEqn,ComputeTimeDerivative)
            
            ComputeA = .FALSE.
         end if
      else
         call this % ComputeJacobian(this % A,dt,time,nEqn,nGradEqn,ComputeTimeDerivative)
         
      end if
      
      call this % A % GetCSR (A)  ! debug
      call A % Visualize('AnJac_Condensed.dat') ! debug
      
      stop ! debug
      
      ! TODO: solve the matrix using SCS
   end subroutine SCS_solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
end module StaticCondensationSolverClass
