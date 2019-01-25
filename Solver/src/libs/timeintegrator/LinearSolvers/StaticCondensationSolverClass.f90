!
!//////////////////////////////////////////////////////
!
!   @File:    StaticCondensationSolverClass.f90
!   @Author:  Andrés Rueda (am.rueda@upm.es)
!   @Created: Tue Dec  4 16:26:02 2018
!   @Last revision date: Fri Jan 25 17:23:15 2019
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 508b6d7bfca8c842ac2d4bdb38ff238e427d2f5c
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
   use DGSEMClass                   , only: DGSem, ComputeTimeDerivative_f
   use GenericLinSolverClass
   use MKLPardisoSolverClass        , only: MKLPardisoSolver_t
   use StaticCondensedMatrixClass   , only: StaticCondensedMatrix_t, INNER_DOF, BOUNDARY_DOF
   use DenseBlockDiagonalMatrixClass, only: DenseBlockDiagMatrix_t
   use CSRMatrixClass               , only: csrMat_t, CSR_MatMatMul, CSR_MatAdd, CSR_MatVecMul
   use Utilities                    , only: AlmostEqual
   implicit none
   
   private
   public StaticCondSolver_t
   
   type, extends(GenericLinSolver_t)   :: StaticCondSolver_t
      type(StaticCondensedMatrix_t)    :: A              ! System matrix
      type(csrMat_t)                   :: Mii_inv        ! Inverse of the inner blocks in CSR (only needed for direct solve)
      type(MKLPardisoSolver_t)         :: linsolver      ! Solver for the condensed system
      real(kind=RP)                    :: Ashift = 0._RP ! Current shift of the Jacobian matrix
      real(kind=RP), allocatable       :: x(:)           ! Solution vector
      real(kind=RP), allocatable       :: bi(:)          ! Right hand side (inner DOFs)
      real(kind=RP), allocatable       :: bb(:)          ! Right hand side ("boundary" DOFs)
   contains
      procedure :: construct           => SCS_construct
      procedure :: destroy             => SCS_destruct
      procedure :: SetOperatorDt       => SCS_SetOperatorDt
      procedure :: ReSetOperatorDt     => SCS_ReSetOperatorDt
      procedure :: solve               => SCS_solve
      procedure :: SetRHS              => SCS_SetRHS
      procedure :: getCondensedSystem  => SCS_getCondensedSystem
      procedure :: getCondensedRHS     => SCS_getCondensedRHS
      procedure :: getGlobalArray      => SCS_getGlobalArray
      procedure :: getLocalArrays      => SCS_getLocalArrays
      procedure :: getSolution         => SCS_getSolution
      procedure :: getX                => SCS_GetX
      procedure :: GetXnorm            => SCS_GetXnorm
      procedure :: GetRnorm            => SCS_GetRnorm
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
      
!     Construct Mii_inv
!     -----------------
      call this % Mii_inv % construct (num_of_Rows = size_i)
      
!
!     Condensed system constructs
!     ---------------------------
      
      call this % linsolver % construct (DimPrb = DimPrb - size_i, MatrixShiftFunc = MatrixShiftFunc, controlVariables = controlVariables)
      
      allocate ( this % bi(size_i) ) 
      allocate ( this % bb(DimPrb - size_i) ) 
      
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
      deallocate (this % bi)
      deallocate (this % bb)
      
      call this % linsolver % destroy
      call this % Mii_inv   % destruct
      
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
      
      if ( AlmostEqual(shift,this % Ashift) ) return
      
      call this % A % Shift(-this % Ashift)
      call this % A % Shift(shift)
      this % Ashift = shift
      
      call this % getCondensedSystem
      
   end subroutine SCS_ReSetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SCS_SetRHS(this, RHS)
      implicit none
      !-arguments-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: RHS(this % DimPrb)
      !-local-variables-----------------------------------------------------
      integer :: i
      !---------------------------------------------------------------------
      
      call this % getLocalArrays(RHS, this % bi, this % bb)
      
   end subroutine SCS_SetRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function SCS_GetX(this) result(x)
      implicit none
      !-arguments-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout) :: this
      real(kind=RP)                            :: x(this % DimPrb)
      !---------------------------------------------------------------------
      
      x = this % x
      
   end function SCS_GetX
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function SCS_GetXnorm(this,TypeOfNorm) result(xnorm)
      implicit none
      !-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout) :: this
      character(len=*)                         :: TypeOfNorm
      real(kind=RP)                            :: xnorm
      !-----------------------------------------------------------
      
      select case (TypeOfNorm)
         case ('infinity')
            xnorm = maxval(abs(this % x))
         case ('l2')
            xnorm = norm2(this % x)
         case default
            stop 'StaticCondensationSolverClass ERROR: Norm not implemented yet'
      end select 
   end function SCS_GetXnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function SCS_GetRnorm(this) result(rnorm)
      implicit none
!
!     ----------------------------------------
!     Currently implemented with infinity norm
!     ----------------------------------------
!
      !-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout) :: this
      real(kind=RP)                            :: rnorm
      !-----------------------------------------------------------
      real(kind=RP)                            :: residual(this % DimPrb)
      !-----------------------------------------------------------
      
!~      residual = this % b - CSR_MatVecMul(this % A, this % x)
!~      rnorm = MAXVAL(ABS(residual))
      !!! TODO: OImplement this!!!
      rnorm = 0._RP
   end function SCS_GetRnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------
!  SCS_getGlobalArray:
!  Get global array from boundary and interior arrays
!  --------------------------------------------------
   subroutine SCS_getGlobalArray(this,x,xi,xb)
      implicit none
      !-arguments-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout)  :: this
      real(kind=RP)            , intent(out)    :: x(this % DimPrb)
      real(kind=RP)            , intent(in)     :: xi(this % A % size_i)
      real(kind=RP)            , intent(in)     :: xb(this % A % size_b)
      !-local-variables-----------------------------------------------------
      integer :: eID    ! Element (block) index
      integer :: i_ind  ! xi index
      integer :: b_ind  ! xb index
      integer :: x_ind  ! x  index
      integer :: i      ! Element DOF index (eq included)
      !---------------------------------------------------------------------
      
      i_ind = 1
      b_ind = 1
      x_ind = 1
      
      do eID = 1, this % A % num_of_Blocks
         do i = 1, this % A % BlockSizes(eID)
            select case (this % A % ElemInfo(eID) % dof_association(i))
               case (INNER_DOF)
                  x(x_ind) = xi(i_ind)
                  i_ind = i_ind + 1
               case (BOUNDARY_DOF)
                  x(x_ind) = xb(b_ind)
                  b_ind = b_ind + 1
            end select
            x_ind = x_ind + 1
         end do  
      end do
      
   end subroutine SCS_getGlobalArray
   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------
!  SCS_getGlobalArray:
!  Get boundary and interior arrays from global array 
!  --------------------------------------------------
   subroutine SCS_getLocalArrays(this,x,xi,xb)
      implicit none
      !-arguments-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: x(this % DimPrb)
      real(kind=RP)            , intent(out)   :: xi(this % A % size_i)
      real(kind=RP)            , intent(out)   :: xb(this % A % size_b)
      !-local-variables-----------------------------------------------------
      integer :: eID    ! Element (block) index
      integer :: i_ind  ! xi index
      integer :: b_ind  ! xb index
      integer :: x_ind  ! x  index
      integer :: i      ! Element DOF index (eq included)
      !---------------------------------------------------------------------
      
      i_ind = 1
      b_ind = 1
      x_ind = 1
      
      do eID = 1, this % A % num_of_Blocks
         do i = 1, this % A % BlockSizes(eID)
            select case (this % A % ElemInfo(eID) % dof_association(i))
               case (INNER_DOF)
                  xi(i_ind) = x(x_ind)
                  i_ind = i_ind + 1
               case (BOUNDARY_DOF)
                  xb(b_ind) = x(x_ind)
                  b_ind = b_ind + 1
            end select
            x_ind = x_ind + 1
         end do  
      end do
      
   end subroutine SCS_getLocalArrays
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
      logical        :: subCompA
      logical        :: mustComputeRHS
      real(kind=RP)  :: xb(this % A % size_b)
      !---------------------------------------------------------------------
      
      mustComputeRHS = .TRUE.
      
!     Compute Jacobian matrix if needed
!     ---------------------------------    
      
      if ( present(ComputeA)) then
         if (ComputeA) then
            call this % ComputeJacobian(this % A,time,nEqn,nGradEqn,ComputeTimeDerivative)
            call this % SetOperatorDt(dt)
            
            ComputeA = .FALSE.
            call this % getCondensedSystem
            mustComputeRHS = .FALSE.
         end if
      else
         call this % ComputeJacobian(this % A,time,nEqn,nGradEqn,ComputeTimeDerivative)
         call this % SetOperatorDt(dt)
         
         call this % getCondensedSystem
         mustComputeRHS = .FALSE.
      end if
      
      if (mustComputeRHS) then 
         call this % getCondensedRHS
      end if
      
      subCompA = .FALSE.
      call this % linsolver % solve( nEqn=nEqn, nGradEqn=nGradEqn, tol = 1.e-6_RP, maxiter=500, time= time, dt=dt, &
                              ComputeTimeDerivative = ComputeTimeDerivative, computeA = subCompA)
      
      xb = this % linsolver % getX()
      call this % getSolution(xb)
      
   end subroutine SCS_solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------
!  SCS_getCondensedSystem:
!  Get condensed system matrix (explicitely: inverting the Mii blocks) and RHS
!  -------------------------------------------------
   subroutine SCS_getCondensedSystem(this)
      implicit none
      !-arguments-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout) :: this
      !-local-variables-----------------------------------------------------
      type(DenseBlockDiagMatrix_t) :: Mii_inv
      type(csrMat_t) :: Mii_inv_Mbi 
      real(kind=RP)  :: Minv_bi(this % A % size_i)
      !---------------------------------------------------------------------
      
      ! Construct auxiliar CSR matrices
      call Mii_inv % construct (num_of_Blocks = this % A % num_of_Blocks)
      call Mii_inv % PreAllocate(nnzs = this % A % inner_blockSizes)
      
      call Mii_inv_Mbi % construct (num_of_Rows = this % A % size_i, num_of_Cols = this % A % size_b)
      
      ! Invert blocks and get CSR
      call this % A % Mii % InvertBlocks_LU (Mii_inv)
      call Mii_inv % getCSR (this % Mii_inv)
      
      ! Get condensed RHS
      Minv_bi = CSR_MatVecMul (this % Mii_inv, this % bi)
      this % linsolver % b = CSR_MatVecMul (this % A % Mib, Minv_bi)
      this % linsolver % b = this % bb - this % linsolver % b
      
      ! Get contensed matrix
      Mii_inv_Mbi = CSR_MatMatMul (this % Mii_inv  , this % A % Mbi)
      this % linsolver % A = CSR_MatMatMul (this % A % Mib, Mii_inv_Mbi )
      this % linsolver % A = CSR_MatAdd (this % A % Mbb, this % linsolver % A, -1._RP)
      
      call this % linsolver % FactorizeJacobian
      
      call Mii_inv % destruct
      call Mii_inv_Mbi % destruct
   end subroutine SCS_getCondensedSystem
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------------
!  SCS_getCondensedRHS:
!  Get condensed RHS (Mii blocks are already inverted)
!  --------------------------------------------------------
   subroutine SCS_getCondensedRHS(this)
      implicit none
      !-arguments-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout) :: this
      !-local-variables-----------------------------------------------------
      type(csrMat_t) :: Mii_inv_Mbi 
      real(kind=RP)  :: Minv_bi(this % A % size_i)
      !---------------------------------------------------------------------
      
      call Mii_inv_Mbi % construct (num_of_Rows = this % A % size_i, num_of_Cols = this % A % size_b)
      
      ! Get condensed RHS
      Minv_bi = CSR_MatVecMul (this % Mii_inv, this % bi)
      this % linsolver % b = CSR_MatVecMul (this % A % Mib, Minv_bi)
      this % linsolver % b = this % bb - this % linsolver % b
      
      call Mii_inv_Mbi % destruct
      
   end subroutine SCS_getCondensedRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------
!  SCS_getSolution:
!  Get global solution out of the boundary solution (xb)
!  -----------------------------------------------------
   subroutine SCS_getSolution(this, xb)
      implicit none
      !-arguments-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: xb(this % A % size_b)
      !-local-variables-----------------------------------------------------
      real(kind=RP) :: xi(this % A % size_i)
      !---------------------------------------------------------------------
      
      xi = CSR_MatVecMul(this % A % Mbi, xb)
      xi = CSR_MatVecMul(this % Mii_inv, this % bi - xi)
      
      call this % getGlobalArray(this % x, xi, xb)
      
   end subroutine SCS_getSolution
end module StaticCondensationSolverClass
