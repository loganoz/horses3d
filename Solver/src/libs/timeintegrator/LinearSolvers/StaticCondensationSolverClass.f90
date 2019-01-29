!
!//////////////////////////////////////////////////////
!
!   @File:    StaticCondensationSolverClass.f90
!   @Author:  Andrés Rueda (am.rueda@upm.es)
!   @Created: Tue Dec  4 16:26:02 2018
!   @Last revision date: Tue Jan 29 18:48:27 2019
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 0f32bff29d29f9d71830bf5971f5e3b189a1d8b8
!
!//////////////////////////////////////////////////////
!
!  StaticCondensationSolverClass: (UNDER DEVELOPMENT!!)
!     Routines for solving Gauss-Lobatto DGSEM representations using static-condensation/substructuring
!
!     -> Only valid for hp-conforming representations (a middle ground can be found -partial condensation)
!
!  TO CHANGE:
!     -> nEqn hard-coded to 5
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
   use StopwatchClass               , only: Stopwatch
   use NodalStorageClass            , only: GAUSSLOBATTO
   use MatrixFreeGMRESClass         , only: MatFreeGMRES_t
   implicit none
   
   private
   public StaticCondSolver_t
   
!
!  ********************************************
!  Main type for the static condensation solver
!  ********************************************
   type, extends(GenericLinSolver_t)   :: StaticCondSolver_t
      type(StaticCondensedMatrix_t)    :: A              ! System matrix
      integer                          :: linsolver      ! Currently used linear solver
      
!     Variables for the direct (sparse LU) solver
!     ------------------------------------------- 
      type(csrMat_t)                   :: Mii_inv           ! Inverse of the inner blocks in CSR (only needed for direct solve)
      type(MKLPardisoSolver_t)         :: pardisoSolver     ! Solver for the condensed system
      
!     Variables for the iterative (GMRES) solver
!     ------------------------------------------
      type(DenseBlockDiagMatrix_t)     :: Mii_LU
      type(MatFreeGMRES_t)             :: gmresSolver
      
      real(kind=RP)                    :: Ashift = 0._RP    ! Current shift of the Jacobian matrix
      real(kind=RP), allocatable       :: x(:)              ! Solution vector
      real(kind=RP), allocatable       :: bi(:)             ! Right hand side (inner DOFs)
      real(kind=RP), allocatable       :: bb(:)             ! Right hand side ("boundary" DOFs)
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
      procedure :: MatrixAction        => SCS_MatrixAction
   end type StaticCondSolver_t
   
!
!  *****************
!  Module parameters
!  *****************
   integer, parameter :: nEqn = 5 ! hard-coded
   
   integer, parameter :: PARDISO = 0
   integer, parameter :: GMRES   = 1
   
!
!  ****************
!  Module variables
!  ****************
   type(StaticCondSolver_t), pointer :: Current_Solver => null()
   
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
      integer :: nelem
      !---------------------------------------------------------------------
!
!     **********************
!     Check needed arguments
!     **********************
!
      if (.not. present(controlVariables)) ERROR stop 'StaticCondSolver_t needs controlVariables'
      if (.not. present(sem)) ERROR stop 'StaticCondSolver_t needs DGSem'
      
      ! TODO: Add conformity check!
!
!     Gauss-Lobatto check
!     -------------------
      
      if (sem % mesh % nodeType /= GAUSSLOBATTO) then
         ERROR stop 'Static Condensation only valid for Gauss-Lobatto discretizations'
      end if
      
      this % DimPrb = DimPrb
      this % p_sem => sem
      
      if ( controlVariables % containsKey("jacobian flag") ) then
         this % JacobianComputation = controlVariables % integerValueForKey("jacobian flag")
      end if
      
      MatrixShift => MatrixShiftFunc
      
!
!     ***************************
!     Construct the system matrix
!     ***************************
!
      nelem  = sem % mesh % no_of_elements
      
      call this % A % construct  (num_of_Rows = DimPrb, &
                                  num_of_Blocks = nelem )
      
!
!     Construct the permutation
!     -------------------------
!
      call this % A % constructPermutationArrays(sem % mesh, nEqn) !,.TRUE. ) ! For ignoring (physical) boundary DOFs 

!
!     ***********
!     Allocations
!     ***********
!
      allocate(this % x(DimPrb))
      allocate ( this % bi(this % A % size_i) ) 
      allocate ( this % bb(DimPrb - this % A % size_i) ) 

!
!     **********************************
!     Construct solver-related variables
!     **********************************
!
      select case ( trim( controlVariables % StringValueForKey("static condensed subsolver",LINE_LENGTH) ) )
         case('pardiso'); this % linsolver = PARDISO
         case('gmres')  ; this % linsolver = GMRES
         case default   ; this % linsolver = PARDISO
      end select
      
      select case (this % linsolver)
         case(PARDISO)
!
!           MKL-Pardiso solver
!           ******************
!
!           Construct Mii_inv
!           -----------------
            call this % Mii_inv % construct (num_of_Rows = this % A % size_i)
      
!
!           Solver
!           ------
            call this % pardisoSolver % construct (DimPrb = DimPrb - this % A % size_i, MatrixShiftFunc = MatrixShiftFunc, controlVariables = controlVariables)
      
         case(GMRES)
!
!           GMRES solver
!           ************
            
            ! Construct auxiliar CSR matrices
            call this % Mii_LU % construct (num_of_Blocks = this % A % num_of_Blocks)
            call this % Mii_LU % PreAllocate(nnzs = this % A % inner_blockSizes)
            
!
!           Solver
!           ------
            call this % gmresSolver % construct  (DimPrb = DimPrb - this % A % size_i, MatrixShiftFunc = MatrixShiftFunc, controlVariables = controlVariables)
            call this % gmresSolver % SetMatrixAction (MatrixAction)
            
      end select
      
      call Stopwatch % CreateNewEvent ("System condensation")
      
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
      
      select case (this % linsolver)
         case(PARDISO)
            call this % pardisoSolver % destroy
            call this % Mii_inv   % destruct
         case(GMRES)
            call this % gmresSolver % destroy
            call this % Mii_LU   % destruct
      end select
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
      class(StaticCondSolver_t), target, intent(inout) :: this
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
      real(kind=RP)  :: xb(this % A % size_b)
      !---------------------------------------------------------------------
      
      Current_Solver => this
      
!     Compute Jacobian matrix if needed
!     ---------------------------------    
      
      if ( present(ComputeA)) then
         if (ComputeA) then
            call this % ComputeJacobian(this % A,time,nEqn,nGradEqn,ComputeTimeDerivative)
            call this % SetOperatorDt(dt)
            
            ComputeA = .FALSE.
            call this % getCondensedSystem
         end if
      else
         call this % ComputeJacobian(this % A,time,nEqn,nGradEqn,ComputeTimeDerivative)
         call this % SetOperatorDt(dt)
         
         call this % getCondensedSystem
      end if
      
      call this % getCondensedRHS
      
      subCompA = .FALSE.
      select case (this % linsolver)
         case (PARDISO)
            call this % pardisoSolver % solve( nEqn=nEqn, nGradEqn=nGradEqn, tol = tol, maxiter=maxiter, time = time, dt=dt, &
                                    ComputeTimeDerivative = ComputeTimeDerivative, computeA = subCompA)
            xb = this % pardisoSolver % getX()
         case (GMRES)
            call this % gmresSolver   % solve( nEqn=nEqn, nGradEqn=nGradEqn, tol = tol, maxiter=maxiter, time = time, dt=dt, &
                                    ComputeTimeDerivative = ComputeTimeDerivative, computeA = subCompA)
            this % niter = this % gmresSolver % niter
            xb = this % gmresSolver % getX()
      end select
      
      call this % getSolution(xb)
      
   end subroutine SCS_solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------
!  SCS_getCondensedSystem:
!  Get condensed system matrix
!  -------------------------------------------------
   subroutine SCS_getCondensedSystem(this)
      implicit none
      !-arguments-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout) :: this
      !-local-variables-----------------------------------------------------
      type(DenseBlockDiagMatrix_t) :: Mii_inv
      type(csrMat_t) :: Mii_inv_Mbi 
      !---------------------------------------------------------------------
      
      select case (this % linsolver)
         case (PARDISO)
            call Stopwatch % Start("System condensation")
            
            ! Construct auxiliar CSR matrices
            call Mii_inv % construct (num_of_Blocks = this % A % num_of_Blocks)
            call Mii_inv % PreAllocate(nnzs = this % A % inner_blockSizes)
            
            call Mii_inv_Mbi % construct (num_of_Rows = this % A % size_i, num_of_Cols = this % A % size_b)
            
            ! Invert blocks and get CSR
            call this % A % Mii % InvertBlocks_LU (Mii_inv)
            call Mii_inv % getTransCSR (this % Mii_inv)
            
            ! Get contensed matrix
            Mii_inv_Mbi = CSR_MatMatMul (this % Mii_inv  , this % A % Mbi, .TRUE.)
            this % pardisoSolver % A = CSR_MatMatMul (this % A % Mib, Mii_inv_Mbi )
            this % pardisoSolver % A = CSR_MatAdd (this % A % Mbb, this % pardisoSolver % A, -1._RP)
            
            call Stopwatch % Pause("System condensation")
            
            call this % pardisoSolver % FactorizeJacobian
            
            call Mii_inv % destruct
            call Mii_inv_Mbi % destruct
         case (GMRES)
            call Stopwatch % Start("System condensation")
            call this % A % Mii % FactorizeBlocks_LU (this % Mii_LU)
            call Stopwatch % Pause("System condensation")
      end select
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
      real(kind=RP)  :: Minv_bi(this % A % size_i)
      !---------------------------------------------------------------------
      
      call Stopwatch % Start("System condensation")
      
      select case (this % linsolver)
         case (PARDISO)
            Minv_bi = CSR_MatVecMul (this % Mii_inv, this % bi, .TRUE.)
            this % pardisoSolver % b = CSR_MatVecMul (this % A % Mib, Minv_bi)
            this % pardisoSolver % b = this % bb - this % pardisoSolver % b
         case (GMRES)
            call this % Mii_LU % SolveBlocks_LU (Minv_bi, this % bi)
            this % gmresSolver % RHS = CSR_MatVecMul (this % A % Mib, Minv_bi)
            this % gmresSolver % RHS = this % bb - this % gmresSolver % RHS
      end select
      
      call Stopwatch % Pause("System condensation")
      
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
      real(kind=RP) :: xi (this % A % size_i)
      !---------------------------------------------------------------------
      call Stopwatch % Pause("Start condensation")
      
      xi = CSR_MatVecMul(this % A % Mbi, xb)
      
      select case (this % linsolver)
         case (PARDISO)
            xi = CSR_MatVecMul(this % Mii_inv, this % bi - xi, .TRUE.)
         case(GMRES)
            call this % Mii_LU % SolveBlocks_LU (xi, this % bi - xi)
      end select
      
      call this % getGlobalArray(this % x, xi, xb)
      
      call Stopwatch % Pause("System condensation")
   end subroutine SCS_getSolution   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------------------
!  SCS_MatrixAction:
!  Perform the product v = Ax, where A is the condensed system matrix
!  --------------------------------------------------------------
   function SCS_MatrixAction(this,x) result(v)
      implicit none
      !-arguments-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: x(this % A % size_b)
      real(kind=RP)                            :: v(this % A % size_b)
      !-local-variables-----------------------------------------------------
      real(kind=RP) :: Mbi_xb (this % A % size_i)   ! Auxiliar vector with size of Mii
      real(kind=RP) :: vi (this % A % size_i)
      real(kind=RP) :: vb (this % A % size_b)
      !---------------------------------------------------------------------
      
      call Stopwatch % Start("System condensation")
      
      ! Compute M_{bi} x_b
      Mbi_xb = CSR_MatVecMul (this % A % Mbi, x)
      
      ! Compute M_{ii}^{-1} M_{bi} x_b
      call this % Mii_LU % SolveBlocks_LU (vi, Mbi_xb)
      
      ! Compute M_{ib} M_{ii}^{-1} M_{bi} x_b
      v = CSR_MatVecMul (this % A % Mib, vi)
      
      ! Compute M_{bb} x_ b
      vb = CSR_MatVecMul (this % A % Mbb, x)
      
      ! Compute M_{bb} x_ b - M_{ib} M_{ii}^{-1} M_{bi} x_b
      v = vb - v
      
      call Stopwatch % Pause("System condensation")
      
   end function SCS_MatrixAction
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  Auxiliar subroutines
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function MatrixAction(x) result(v)
      implicit none
      !-arguments-----------------------------------------------------------
      real(kind=RP)            , intent(in)    :: x(:)
      real(kind=RP)                            :: v(size(x))
      !---------------------------------------------------------------------
      
      v = Current_Solver % MatrixAction(x)
      
   end function MatrixAction
end module StaticCondensationSolverClass
