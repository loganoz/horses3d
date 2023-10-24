!//////////////////////////////////////////////////////
!
!  StaticCondensationSolverClass:
!     Routines for solving Gauss-Lobatto DGSEM representations using static-condensation/substructuring
!
module StaticCondensationSolverClass
   use SMConstants
   use DGSEMClass                   , only: DGSem, ComputeTimeDerivative_f
   ! Linear solvers
   use GenericLinSolverClass
   use MKLPardisoSolverClass        , only: MKLPardisoSolver_t
   use PetscSolverClass             , only: PetscKspLinearSolver_t
   ! Matrices
   use GenericMatrixClass           , only: Matrix_t
   use StaticCondensedMatrixClass   , only: StaticCondensedMatrix_t, INNER_DOF, BOUNDARY_DOF, SC_MATRIX_CSR, SC_MATRIX_PETSC
   use DenseBlockDiagonalMatrixClass, only: DenseBlockDiagMatrix_t
   use CSRMatrixClass               , only: csrMat_t, CSR_MatAdd
   use PETScMatrixClass             , only: PETSCMatrix_t
   ! Extras
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
      
!     Variables for the matrix (pardiso/PETSc) solvers
!     ------------------------------------------------
      class(Matrix_t)          , allocatable :: Mii_inv           ! Inverse of the inner blocks in CSR (only needed for matrix solvers)
      class(GenericLinSolver_t), allocatable :: matSolver         ! Solver for the condensed system
      
!     Variables for the matrix-free (GMRES) solver
!     --------------------------------------------
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
   
!
!  Sub solvers
!  -----------
   integer, parameter :: SSOLVER_PARDISO    = 0
   integer, parameter :: SSOLVER_PETSC      = 1
   integer, parameter :: SSOLVER_MATF_GMRES = 2
   
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
   subroutine SCS_construct(this,DimPrb, globalDimPrb, nEqn,controlVariables,sem,MatrixShiftFunc)
      implicit none
      !-arguments-----------------------------------------------------------
      class(StaticCondSolver_t), intent(inout), target :: this
      integer                  , intent(in)            :: DimPrb
      integer                  , intent(in)            :: globalDimPrb
      integer                  , intent(in)            :: nEqn
      type(FTValueDictionary)  , intent(in), optional  :: controlVariables
      type(DGSem), target                  , optional  :: sem
      procedure(MatrixShift_FCN)                       :: MatrixShiftFunc
      !-local-variables-----------------------------------------------------
      integer :: nelem
      character(len=LINE_LENGTH) :: linsolver
      integer :: MatrixType
      !---------------------------------------------------------------------
      
      call this % GenericLinSolver_t % construct(DimPrb,globalDimPrb, nEqn,controlVariables,sem,MatrixShiftFunc)
!
!     **********************
!     Check needed arguments
!     **********************
!
      if (.not. present(controlVariables)) error stop 'StaticCondSolver_t needs controlVariables'
      if (.not. present(sem)) error stop 'StaticCondSolver_t needs DGSem'
      
      ! TODO: Add conformity check?
!
!     Gauss-Lobatto check
!     -------------------
      
      if (sem % mesh % nodeType /= GAUSSLOBATTO) then
         error stop 'Static Condensation only valid for Gauss-Lobatto discretizations'
      end if
      
      this % DimPrb = DimPrb
      this % p_sem => sem
       
      MatrixShift => MatrixShiftFunc
!
!     ****************
!     Select subsolver
!     ****************
!
      linsolver = trim( controlVariables % StringValueForKey("static condensed subsolver",LINE_LENGTH) )
      if ( trim(linsolver) == '' ) linsolver = 'pardiso' ! default value
      
      select case ( trim(linsolver) )
         case('pardiso')
            this % linsolver = SSOLVER_PARDISO
            allocate(csrMat_t :: this % Mii_inv)
            allocate(MKLPardisoSolver_t :: this % matSolver)
            MatrixType = SC_MATRIX_CSR
         case('petsc')
            this % linsolver = SSOLVER_PETSC
            allocate(PETSCMatrix_t :: this % Mii_inv)
            allocate(PetscKspLinearSolver_t :: this % matSolver)
            MatrixType = SC_MATRIX_PETSC
         case('gmres')
            this % linsolver = SSOLVER_MATF_GMRES
         case default   
            write(*,'(A)') 'Not recognized static condensed subsolver.'
            write(*,'(A)') 'Options are:'
            write(*,'(A)') '  * pardiso'
            write(*,'(A)') '  * petsc'
            write(*,'(A)') '  * gmres (matrix-free)'
            error stop
      end select
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
      call this % A % constructPermutationArrays(sem % mesh, nEqn, MatrixType) !,.TRUE. ) ! For ignoring (physical) boundary DOFs 

!
!     ***********
!     Allocations
!     ***********
!
      allocate(this % x(DimPrb))
      allocate ( this % bi(this % A % size_i) ) 
      allocate ( this % bb(DimPrb - this % A % size_i) ) 
!
!     ***********************
!     Construct linear solver
!     ***********************
!     
      select case (this % linsolver)
         case(SSOLVER_PARDISO, SSOLVER_PETSC)
!
!           Matrix solvers
!           **************

            ! Construct solver
            call this % matSolver % construct(DimPrb = DimPrb - this % A % size_i, &
                                              globalDimPrb = globalDimPrb - this % A % size_i, & ! change this for MPI
                                              nEqn = nEqn, &
                                              MatrixShiftFunc = MatrixShiftFunc, &
                                              controlVariables = controlVariables)
            
            ! This is for adding the block sizes to the subsolver Jacobian definition (only valid in sequential)
            allocate ( this % matSolver % Jacobian % ndofelm_l(nelem) )
            this % matSolver % Jacobian % ndofelm_l = this % A % BlockSizes - this % A % inner_blockSizes
            
            ! Construct auxiliary matrix (nor really needed now since the matrix is constructed in this % A % getSchurComplement())
!~            call this % Mii_inv % construct (num_of_Rows = this % A % size_i)
            
         case(SSOLVER_MATF_GMRES)
!
!           Matrix-free solver
!           ******************
            
            ! Construct auxiliary CSR matrices
            call this % Mii_LU % construct (num_of_Blocks = this % A % num_of_Blocks)
            call this % Mii_LU % PreAllocate(nnzs = this % A % inner_blockSizes)
            
            ! Construct solver
            call this % gmresSolver % construct (DimPrb = DimPrb - this % A % size_i, &
                                                 globalDimPrb = globalDimPrb - this % A % size_i, & ! change this for MPI
                                                 nEqn = nEqn, &
                                                 MatrixShiftFunc = MatrixShiftFunc, &
                                                 controlVariables = controlVariables)
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
!
!        Matrix solvers
!        **************
         case(SSOLVER_PARDISO, SSOLVER_PETSC)
!
!           Destroy auxiliary matrix (must be done prior to destroying solver)
!           -----------------------------------------------------------------
            call this % Mii_inv   % destruct
            deallocate (this % Mii_inv)
!
!           Destroy solver
!           --------------
            call this % matSolver % destroy
            deallocate (this % matSolver)
!
!        Matrix-free solvers
!        *******************
         case(SSOLVER_MATF_GMRES)
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
            error stop 'StaticCondensationSolverClass ERROR: Norm not implemented yet'
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
      
      select case (this % linsolver)
         case (SSOLVER_MATF_GMRES)
            rnorm = this % gmresSolver % GetRnorm()
         case default
            rnorm = this % matSolver % GetRnorm()
      end select
      
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
            call this % Jacobian % Compute (this % p_sem, nEqn, time, this % A, ComputeTimeDerivative)
            call this % SetOperatorDt(dt)
            
            ComputeA = .FALSE.
            call this % getCondensedSystem
         end if
      else
         call this % Jacobian % Compute (this % p_sem, nEqn, time, this % A, ComputeTimeDerivative)
         call this % SetOperatorDt(dt)
         
         call this % getCondensedSystem
      end if
      
      call this % getCondensedRHS
      
      subCompA = .FALSE.
      
      select case (this % linsolver)
         case (SSOLVER_PARDISO, SSOLVER_PETSC)
            call this % matSolver % solve( nEqn=nEqn, nGradEqn=nGradEqn, tol = tol, maxiter=maxiter, time = time, dt=dt, &
                                    ComputeTimeDerivative = ComputeTimeDerivative, computeA = subCompA)
            this % niter = this % matSolver % niter
            xb = this % matSolver % getX()
         case (SSOLVER_MATF_GMRES)
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
      !---------------------------------------------------------------------
      
      select case (this % linsolver)
         case (SSOLVER_PARDISO, SSOLVER_PETSC)
            call Stopwatch % Start("System condensation")
            
            ! Construct auxiliary matrix for Mii⁻¹
            call Mii_inv % construct (num_of_Blocks = this % A % num_of_Blocks)
            call Mii_inv % PreAllocate(nnzs = this % A % inner_blockSizes)
            
            ! Invert blocks and get corresponding sparse matrix
            call this % A % Mii % InvertBlocks_LU (Mii_inv)
            call this % Mii_inv % ConstructFromDiagBlocks(Mii_inv % num_of_Blocks, Mii_inv % Blocks, Mii_inv % BlockIdx, Mii_inv % BlockSizes)
            call Mii_inv % destruct
            
            ! Additional operations (solver specific)
            select type (matSolver => this % matSolver)
               class is (PetscKspLinearSolver_t)
                  ! Get condensed matrix
                  call this % A % getSchurComplement(this % Mii_inv, matSolver % A)
                  call Stopwatch % Pause("System condensation")
                  
                  ! Set PETSc matrix and preconditioner
                  call matSolver % SetPreconditioner
                  
               class is (MKLPardisoSolver_t)
                  ! Get condensed matrix
                  call this % A % getSchurComplement(this % Mii_inv, matSolver % A)
                  call Stopwatch % Pause("System condensation")
                  
                  ! Factorize Jacobian
                  call matSolver % FactorizeJacobian
            end select
            
         case (SSOLVER_MATF_GMRES)
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
      real(kind=RP)  :: bb     (this % A % size_b)
      !---------------------------------------------------------------------
      
      call Stopwatch % Start("System condensation")
      
      select case (this % linsolver)
         case (SSOLVER_PARDISO, SSOLVER_PETSC)
            Minv_bi = this % Mii_inv % MatVecMul (this % bi)
            bb = this % A % Mib % MatVecMul (Minv_bi)
            call this % matSolver % SetRHS (this % bb - bb)
         case (SSOLVER_MATF_GMRES)
            call this % Mii_LU % SolveBlocks_LU (Minv_bi, this % bi)
            this % gmresSolver % RHS = this % A % Mib % MatVecMul (Minv_bi)
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
      call Stopwatch % Start("System condensation")
      
      xi = this % A % Mbi % MatVecMul (xb)
      
      select case (this % linsolver)
         case (SSOLVER_PARDISO, SSOLVER_PETSC)
            xi = this % Mii_inv % MatVecMul (this % bi - xi)
         case(SSOLVER_MATF_GMRES)
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
      real(kind=RP) :: Mbi_xb (this % A % size_i)   ! Auxiliary vector with size of Mii
      real(kind=RP) :: vi (this % A % size_i)
      real(kind=RP) :: vb (this % A % size_b)
      !---------------------------------------------------------------------
      
      call Stopwatch % Start("System condensation")
      
      ! Compute M_{bi} x_b
      Mbi_xb = this % A % Mbi % MatVecMul (x)
      
      ! Compute M_{ii}^{-1} M_{bi} x_b
      call this % Mii_LU % SolveBlocks_LU (vi, Mbi_xb)
      
      ! Compute M_{ib} M_{ii}^{-1} M_{bi} x_b
      v = this % A % Mib % MatVecMul (vi)
      
      ! Compute M_{bb} x_ b
      vb = this % A % Mbb % MatVecMul (x)
      
      ! Compute M_{bb} x_ b - M_{ib} M_{ii}^{-1} M_{bi} x_b
      v = vb - v
      
      call Stopwatch % Pause("System condensation")
      
   end function SCS_MatrixAction
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  Auxiliary subroutines
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