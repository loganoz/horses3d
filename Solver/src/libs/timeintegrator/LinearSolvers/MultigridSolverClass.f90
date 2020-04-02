!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!      MultigridSolverClass.f90
!      Created: Wed Nov 6 17:45:26 CET 2019
!      Version: 1.0 (Wed Nov 6 17:45:26 CET 2019)
!      Author: Wojciech Laskowski (wj.laskowski@upm.es) based on an old, unfinished routine by Andrés Rueda (XXXXX - last commit before deleting)
!
!      Class for solving a linear system obtained from a DGSEM discretization using p-Multigrid.
!
!      Control variables:  
!      no_levels :: number of MG levels, IF NOT (no_levels = 2)        
!      define levels x :: hard define N on each level, IF NOT \Delta N_{x} = 1
!      define levels y ...
!      define levels z ...
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------
!  TODO:
!  ----------------
!  1. Deallacotion: 
!     a. DGSem on each level
!     b. Jac on each level
!     c. Prol/Rest operators 
!  2. Generalize MG to operate on nonconforming p-mesh (elements with different pol. orders)
!  3. Add more smoothers - so far only diag-jacobi implemented.
!  ----------------
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
module MultigridSolverClass
   use GenericLinSolverClass
   use CSRMatrixClass
   use SMConstants
   use PhysicsStorage
   use StorageClass
   use PolynomialInterpAndDerivsModule
   use GaussQuadrature
   use DGSEMClass
   use TimeIntegratorDefinitions
   use MatrixClass
   use NumericalJacobian      , only: NumJacobian_t
   use AnalyticalJacobian     , only: AnJacobian_t
   use JacobianComputerClass  , only: JacobianComputer_t, GetJacobianFlag
   use GenericSmoother

   implicit none
   
   private
   public MultigridSolver_t, TemporaryElementStorage_t
!
!  ------------------------------------------------
!  Local temporary element storage.
!  ------------------------------------------------
!
   type :: TemporaryElementStorage_t
!-----Variables-----------------------------------------------------------
      real(kind=RP), dimension(:,:,:,:), allocatable :: X 
      real(kind=RP), dimension(:,:,:,:), allocatable :: R 
   end type TemporaryElementStorage_t
!
!  ------------------------------------------------
!  Multigrid type (Linear Solver class extension)
!  ------------------------------------------------
!
   type, extends(GenericLinSolver_t) :: MultigridSolver_t
!-----Variables-----------------------------------------------------------
      type(csrMat_t)                                :: A                     ! Matrix to solve
      type(DenseBlockDiagMatrix_t),     allocatable :: Atmp                  ! tmp
      type(BJSmooth_t)            ,     allocatable :: BJSmoother            ! smoother
      real(kind=RP), dimension(:) ,     allocatable :: x                     ! Solution vector
      real(kind=RP), dimension(:) ,     allocatable :: b                     ! Right hand side
      real(kind=RP), dimension(:) ,     allocatable :: r                     ! Residual
      real(kind=RP)                                 :: rnorm                 ! L2 norm of residual
      real(kind=RP)                                 :: tol                   ! Tolerance
      real(kind=RP)                                 :: maxiter               ! Max. # iterations
      integer                                       :: MGlevel               ! Current level
      integer                                       :: Nx                    ! Polynomial in X
      integer                                       :: Ny                    ! Polynomial in Y
      integer                                       :: Nz                    ! Polynomial in Z
      type(MultigridSolver_t), pointer              :: Child                 ! Coarser level: MGlevel-1
      type(MultigridSolver_t), pointer              :: Parent                ! Finer level: MGlevel+1
      type(TemporaryElementStorage_t) , allocatable :: LocalStorage(:)       ! Storage
      ! 1D prolongation/restriction operators in each direction 
      real(kind=RP), dimension(:,:) ,   allocatable :: RestX(:,:)             
      real(kind=RP), dimension(:,:) ,   allocatable :: ProlX(:,:)             
      real(kind=RP), dimension(:,:) ,   allocatable :: RestY(:,:)             
      real(kind=RP), dimension(:,:) ,   allocatable :: ProlY(:,:)             
      real(kind=RP), dimension(:,:) ,   allocatable :: RestZ(:,:)             
      real(kind=RP), dimension(:,:) ,   allocatable :: ProlZ(:,:)             
      ! 3D prolongation/restriction operators in each direction 
      real(kind=RP), dimension(:,:) ,   allocatable :: Rest3D(:,:)             
      real(kind=RP), dimension(:,:) ,   allocatable :: Prol3D(:,:)             
      ! test for the Jacobian 
      real(kind=RP), dimension(:,:) ,   allocatable :: RestJacX(:,:)
      real(kind=RP), dimension(:,:) ,   allocatable :: RestJacY(:,:)
      real(kind=RP), dimension(:,:) ,   allocatable :: RestJacZ(:,:)
      
   contains
!-----Subroutines-----------------------------------------------------------
      procedure :: Construct            => MG_Construct 
      procedure :: Destroy              => MG_Destruct 
      procedure :: SetRHS               => MG_SetRHS
      procedure :: GetX                 => MG_GetX
      procedure :: Solve                => MG_Solve 
      procedure :: Getrnorm             => MG_Getrnorm
      procedure :: Getxnorm             => MG_Getxnorm
      procedure :: SetRHSValue          => MG_SetRHSValue 
      procedure :: SetRHSValues         => MG_SetRHSValues 
      procedure :: MG_CreateProlongationOperator
      procedure :: MG_CreateRestrictionOperator
      !procedure :: MG_Create3DProlongationMatrix
      !procedure :: MG_Create3DRestrictionMatrix
      !procedure :: prolong              => MG_Prolongation
      !procedure :: restrict             => MG_Restriction
      procedure :: MG_Prolongation
      procedure :: MG_Restriction

   end type MultigridSolver_t
!
!  ----------------
!  Module variables
!  ----------------
!
   integer                            :: no_levels    ! Total number of multigrid levels
   integer, dimension(:), allocatable :: MG_levels_x  ! multigrid levels
   integer, dimension(:), allocatable :: MG_levels_y  ! multigrid levels
   integer, dimension(:), allocatable :: MG_levels_z  ! multigrid levels
   integer, dimension(:), allocatable :: pre_smooths  ! pre-smoothing operations
   integer, dimension(:), allocatable :: pos_smooths  ! post-smoothing operations
   integer                            :: nelem        ! Number of elements
   integer                            :: S_SMOOTHER   ! Smoother type
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Construct(this,DimPrb,globalDimPrb,nEqn,controlVariables,sem,MatrixShiftFunc)
      implicit none
!-----Arguments-----------------------------------------------------------
      class(MultigridSolver_t),  intent(inout), target :: this
      integer                  , intent(in)            :: DimPrb
      integer                  , intent(in)            :: globalDimPrb        
      integer                  , intent(in)            :: nEqn
      type(FTValueDictionary)  , intent(in), optional  :: controlVariables
      type(DGSem), target                  , optional  :: sem
      procedure(MatrixShift_FCN)                       :: MatrixShiftFunc     ! TODO: Make this optional
      procedure(ComputeTimeDerivative_f)               :: ComputeTimeDerivative
!-----Local-Variables-----------------------------------------------------------
      character(len=LINE_LENGTH)                 :: pc
      integer                                    :: i
      real(kind=RP)                              :: tmp_1

      ! Check dims, MPI and allocate Jacobian
      call this % GenericLinSolver_t % construct(DimPrb, globalDimPrb, nEqn,controlVariables,sem,MatrixShiftFunc)
      MatrixShift => MatrixShiftFunc ! FIXME: ask whether we need this 

!     ***********************************         
!     Set variables from controlVariables
!     ***********************************
      if (.not. present(controlVariables)) stop 'Fatal error: MultigridSolver needs controlVariables.'
      if (present(controlVariables)) then 

        ! Multigrid coarse grids
        ! *********************************************************
        if ( controlVariables % containsKey("multigrid levels") ) then
            no_levels = controlVariables % integerValueForKey("multigrid levels")
        else
            no_levels = 2
        end if

        print *, "No of pMG levels: ", no_levels
        allocate(MG_levels_x(no_levels))
        allocate(MG_levels_y(no_levels))
        allocate(MG_levels_z(no_levels))
        allocate(pre_smooths(no_levels))
        allocate(pos_smooths(no_levels))

        if ( controlVariables % containsKey("multigrid levels") ) then
            ! Levels in X
            pc = controlVariables % StringValueForKey("define levels x",LINE_LENGTH)
            do i = 1, len_trim(pc)
              read(pc(i:i),'(i1)') MG_levels_x(i)
            end do
            ! Levels in Y
            pc = controlVariables % StringValueForKey("define levels y",LINE_LENGTH)
            do i = 1, len_trim(pc)
              read(pc(i:i),'(i1)') MG_levels_y(i)
            end do
            ! Levels in Z
            pc = controlVariables % StringValueForKey("define levels z",LINE_LENGTH)
            do i = 1, len_trim(pc)
              read(pc(i:i),'(i1)') MG_levels_z(i)
            end do
        else
            ERROR stop ':: FIXME: Need to finish this.'
        end if
        ! *********************************************************

        ! Smoother
        ! *********************************************************
        if ( controlVariables % containsKey("smoother") ) then
        select case ( trim( controlVariables % StringValueForKey("smoother",LINE_LENGTH) ) )
           case ('point-jacobi')
             S_SMOOTHER = S_POINTJAC
             print *, 'GenericSmoother :: WARNING p. Jacobi smoother is not suitable for DGSEM Jacobians'
           case ('block-jacobi')
             S_SMOOTHER = S_BLOCKJAC
           case default
             ERROR Stop "MultigridSolver :: Wrong smoother "
        end select
        end if

        if ( controlVariables % containsKey("pre smooths") ) then
          pc = controlVariables % StringValueForKey("pre smooths",LINE_LENGTH)
          pre_smooths = getnosmooths(pc,no_levels)
        print *, "No of pre-smoothing operations: ", pre_smooths
        end if

        if ( controlVariables % containsKey("post smooths") ) then
          pc = controlVariables % StringValueForKey("post smooths",LINE_LENGTH)
          pos_smooths = getnosmooths(pc,no_levels)
        print *, "No of post-smoothing operations: ", pos_smooths
        end if
        ! *********************************************************
      end if 
      nelem = size(sem % mesh % elements)

      this % p_sem => sem ! this is for generic linear solver class
      ! this % sem_lvl => sem ! this is for local sem on each level
      this % DimPrb = DimPrb

      print *, "My fine mesh: "
      print *, " ************************** "
      print *, "Mesh total DOF: ", this % p_sem % mesh % NDOF
      print *, "DimPrb: ", this % DimPrb 
      print *, " ************************** "
  
      ! Construct initial variables for the finest level and call recursive constructor
      ! *********************************************************
      this % MGlevel = no_levels
      this % Nx = MG_levels_x(no_levels)
      this % Ny = MG_levels_y(no_levels)
      this % Nz = MG_levels_z(no_levels)

      ! MAINJAC
      call this % A % construct(num_of_Rows = DimPrb, withMPI = .false.)
      call this % Jacobian % Configure (sem % mesh, nEqn, this % A)
      ! MAINJAC
      

      call MG_Levels_Construct(this,no_levels,controlVariables,nEqn)
      ! *********************************************************


   end subroutine MG_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine MG_Levels_Construct(Me,lvl,controlVariables,nEqn)
      implicit none
!-----Arguments-----------------------------------------------------------
      type(MultigridSolver_t),  intent(inout), target :: Me
      integer,                  intent(in)            :: nEqn
      integer,                  intent(in)            :: lvl
      class(FTValueDictionary), intent(in)            :: controlVariables                
!-----Local-Variables-----------------------------------------------------------
      type(MultigridSolver_t), pointer  :: Child_p          ! Pointer to Child
      integer                           :: i,j,k
      integer                           :: nnzs(nelem)

      ! character(len=*) :: meshFileName_
      !integer, allocatable              :: Nx(:), Ny(:), Nz(:)
      ! integer                           :: polynomialOrder(3)
      logical                           :: success            

      ! for coarse jacobian
      integer                                 :: JacobianComputation ! Get Jacobian type 

      allocate ( Me % b(Me % DimPrb) ) 
      allocate ( Me % x(Me % DimPrb) ) 
      allocate ( Me % r(Me % DimPrb) ) 

      ! Construct smoother
      select case (S_SMOOTHER) 
      case (S_POINTJAC)
         print *, "We're doing PJ"
      case (S_BLOCKJAC)
         ! Allocate and construct Block-Jacobian matrix
         allocate (Me % Atmp)
         call Me % Atmp % construct (num_of_Blocks = Me % p_sem % mesh % no_of_elements)
         do j=1,nelem
             nnzs(j) = nEqn*(Me%Nx+1)*(Me%Ny+1)*(Me%Nz+1) 
         end do
         call Me % Atmp % PreAllocate (nnzs=nnzs)

         ! Allocate smoother
         allocate ( Me % BJSmoother )
         call Me % BJSmoother % Construct ( Me % Atmp )
      case default
         ERROR Stop "Shouldnt be here"
      end select 



      ALLOCATE ( Me % LocalStorage(nelem))
! !$omp parallel do private(Q1,Q2,Q3,Q4) schedule(runtime)
      do k = 1, nelem
         allocate(Me % LocalStorage(k) % X (nEqn,0:Me%Nx,0:Me%Ny,0:Me%Nz))
         allocate(Me % LocalStorage(k) % R (nEqn,0:Me%Nx,0:Me%Ny,0:Me%Nz))
         Me % LocalStorage(k) % X = 0._RP
         Me % LocalStorage(k) % R = 0._RP
      end do   
! !$omp end parallel do

      ! Construct coarser level
      ! *********************************************************
      if (lvl > 1) then
        ! Allocate the child
        ! *********************************************************
        allocate(Me % Child) ! allocate coarser MG level (child)
        Child_p => Me % Child ! set local pointer to a child
        Me % Child % Parent => Me ! set child's parent pointer to this level (Me)
        Child_p % MGlevel = lvl - 1
        Child_p % Nx = MG_levels_x(lvl - 1)
        Child_p % Ny = MG_levels_y(lvl - 1)
        Child_p % Nz = MG_levels_z(lvl - 1)
        ! *********************************************************

        ! Create coarse DGSem
        ! *********************************************************
        ! polynomialOrder = (/ Child_p%Nx, Child_p%Ny, Child_p%Nz /)
        allocate (Child_p % p_sem)
        print *, "Just to make sure... "
        call Child_p % p_sem % construct (  controlVariables  = controlVariables, &
        polynomialOrder = (/ Child_p%Nx, Child_p%Ny, Child_p%Nz /), success           = success)
        if (.not. success) ERROR STOP "Multigrid: Problem creating coarse solver."

        Child_p % DimPrb = Child_p % p_sem % mesh % NDOF * nEqn 
        print *, "My coarse mesh: "
        print *, " ************************** "
        print *, "Mesh total DOF: ", Child_p % p_sem % mesh % NDOF
        print *, "DimPrb: ", Child_p % DimPrb 
        print *, "DOF of the each element: "
        ! do i = 1, size(Child_p % p_sem % mesh % elements) 
        !    print *, i, Child_p % p_sem % mesh % elements(i) % storage % NDOF
        ! end do
        print *, " ************************** "
        ! *********************************************************

        ! Construct Jacobian on coarse level
        ! *********************************************************
        call Child_p % A % construct(num_of_Rows = Child_p % DimPrb, withMPI = .false.)
        JacobianComputation = GetJacobianFlag()
        print *, "Jacobian Type: ", JacobianComputation
        select case (JacobianComputation)
          case (ANALYTICAL_JACOBIAN) ; allocate(AnJacobian_t  :: Child_p % Jacobian)
          case (NUMERICAL_JACOBIAN ) ; allocate(NumJacobian_t :: Child_p % Jacobian)
          case default
             ERROR stop 'Invalid jacobian type. FIXME: '
        end select

        call Child_p % Jacobian % construct(Child_p % p_sem % mesh, nEqn)
        call Child_p % Jacobian % configure(Child_p % p_sem % mesh, nEqn, Child_p % A)
        ! call Child_p % A % Visualize("JacC.txt") ! write Jacobian to a file
        ! *********************************************************

        ! Prolongation restriction operators
        ! *********************************************************
        allocate (Me  % RestX(0:Child_p%Nx,0:Me%Nx)                 )
        allocate (Child_p % ProlX(0:Me%Nx,0:Child_p%Nx)             )
        allocate (Me  % RestY(0:Child_p%Ny,0:Me%Ny)                 )
        allocate (Child_p % ProlY(0:Me%Ny,0:Child_p%Ny)             )
        allocate (Me  % RestZ(0:Child_p%Nz,0:Me%Nz)                 )
        allocate (Child_p % ProlZ(0:Me%Nz,0:Child_p%Nz)             )
        ! Test for the Jacobian
        allocate (Me  % RestJacX(0:Child_p%Nx,0:Me%Nx)                 )
        allocate (Me  % RestJacY(0:Child_p%Nx,0:Me%Nx)                 )
        allocate (Me  % RestJacZ(0:Child_p%Nx,0:Me%Nx)                 )

        ! 3D Prolongation restriction operators
        ! *********************************************************
        allocate (Me      % Rest3D( 1:( (Child_p%Nx+1) * (Child_p%Ny+1) * (Child_p%Nz+1) ), 1:( (Me%Nx+1) * (Me%Ny+1) * (Me%Nz+1) ) )  )
        allocate (Child_p % Prol3D( 1:( (Me%Nx+1) * (Me%Ny+1) * (Me%Nz+1) ), 1:( (Child_p%Nx+1) * (Child_p%Ny+1) * (Child_p%Nz+1) ) )  )
        ! allocate (Child_p % Prol3D(0:Me%Nx,0:Child_p%Nx)  )

        print *, "My interpolation matrices: "
        call MG_CreateRestrictionOperator  ( Me , .true.)
        ! call MG_CreateRestrictionOperator  ( Me , .false.)
        ! call MG_CreateRestrictionOperator  ( Me )
        print *, "Level: ", lvl
        call MG_CreateProlongationOperator ( Child_p )

        ! *********************************************************

        ! recursive call
        call MG_Levels_Construct(Me % Child,lvl-1,controlVariables,nEqn)
      end if
      ! *********************************************************

      ! Printing
      ! *********************************************************
      print *, " ************************** "
      print *, "We're on level: ", Me%MGlevel
      print *, "Pol. orders: in x/y/z:", Me%Nx, "/", Me%Ny, "/", Me%Nz
      print *, " ************************** "
      ! *********************************************************
   end subroutine MG_Levels_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Destruct(this)
     implicit none
     class(MultigridSolver_t), intent(inout) :: this

     call this % A % destruct

     deallocate(MG_levels_x)
     deallocate(MG_levels_y)
     deallocate(MG_levels_z)
   end subroutine MG_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Solve(this, nEqn, nGradEqn, ComputeTimeDerivative, tol, maxiter, time, dt, ComputeA)
      implicit none
      class(MultigridSolver_t), target, intent(inout) :: this
      type(MultigridSolver_t) , pointer               :: pMG     
      integer,       intent(in)                       :: nEqn, nGradEqn
      procedure(ComputeTimeDerivative_f)              :: ComputeTimeDerivative
      real(kind=rp), optional                         :: tol
      integer      , optional                         :: maxiter
      real(kind=rp), optional                         :: time
      real(kind=rp), optional                         :: dt
      logical      , optional , intent(inout)         :: ComputeA
      !---------------------------------------------------------------------
      class(csrMat_t), pointer                        :: pA
      integer                                         :: niter
      real(kind=rp)                                   :: tmpsize
      integer                                         :: i, j, k
      character(len=1024)                             :: filename
      logical                                         :: file_exists

      if ( present(ComputeA)) then
         if (ComputeA) then
            print *, "Computing Jacobian... "
            call MG_ComputeJacobians( this,no_levels,ComputeTimeDerivative,Time,nEqn )
            print *, "   ... done. "
            ComputeA = .FALSE.
         end if
      else
         print *, "Computing Jacobian... "
         call MG_ComputeJacobians( this,no_levels,ComputeTimeDerivative,Time,nEqn )
         print *, "   ... done. "
      end if

      ! ! some printing 
      ! print *, "Elements Q: "
      ! print *, "-----------------------------"
      ! do k = 1, nelem 
      !    print *, k, " out of ", nelem
      !    print *, this % p_sem % mesh % elements(k) % storage % Q
      ! end do
      ! print *, "-----------------------------"
      ! print *, "Global Q: "
      ! print *, "-----------------------------"
      ! print *, this % p_sem % mesh % storage % Q 
      ! print *, "-----------------------------"


      call this % child % p_sem % mesh % storage % local2globalq (this % child % p_sem % mesh % storage % NDOF)
      ! print *, this % p_sem % mesh % storage % NDOF
      ! print *, this % child % p_sem % mesh % storage % NDOF
      open(66, file = 'Q_tot.dat')
      write(66,'(1E19.11)') this % child % p_sem % mesh % storage % Q
      close(66) 

      open(67, file = 'Q_els.dat')
      do k = 1, nelem 
         write(67,*) "Element ", k
         write(67,'(1E19.11)') this % child % p_sem % mesh % elements(k) % storage % Q
      end do
      close(67) 

      print *, "Q1"
      print *, this % child % p_sem % mesh % elements(1) % storage % Q(1,:,:,:)
      print *, "Q2"
      print *, this % child % p_sem % mesh % elements(1) % storage % Q(2,:,:,:)
      print *, "Q3"
      print *, this % child % p_sem % mesh % elements(1) % storage % Q(3,:,:,:)
      print *, "Q4"
      print *, this % child % p_sem % mesh % elements(1) % storage % Q(4,:,:,:)
      print *, "Q5"
      print *, this % child % p_sem % mesh % elements(1) % storage % Q(5,:,:,:)

      print *, size(this%Rest3D,1), size(this%Rest3D,2)
      print *, size(this%child%Prol3D,1), size(this%child%Prol3D,2)
      
      open(69, file = 'Prol3D.dat')
      write(69,'(1E19.11)') this % child % Prol3D 
      close(69) 

      open(70, file = 'Rest3D.dat')
      write(70,'(1E19.11)') this % Rest3D 
      close(70) 

      ! error STOP "TBC Wojtek"

      ! ! Visualize RHS
      ! open(1, file = 'b.txt', status = 'new')  
      ! write(1,*) this % b  
      ! close(1) 
      ! Error Stop "TBC Wojtek"

      this % niter = 0
      this % maxiter = 30
      this % converged = .false.
      this % tol = 1e-6
     
      this % x = 0.0_RP ! setting initial solution
      do while ( (.not. this % converged ) .and. (this % niter .lt. this % maxiter) )
         call MG_VCycle( this, no_levels, nEqn)

         this % r = CSR_MatVecMul( this % A, this % x ) 
         do i = 1 , this % DimPrb
            this % r(i) = this % b(i) - this % r(i)
         end do 
         this % rnorm = norm2(this % r)
         this % niter = this % niter + 1
         print *, "Cycle ", this % niter, ". Res: ", this % rnorm
         if (this % rnorm < this % tol) this % converged = .true.
      end do

      ERROR stop ':: TBC Wojtek'
   end subroutine MG_Solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine MG_VCycle(Me,lvl,nEqn)
      implicit none
      !-----------------------------------------------------
      type(MultigridSolver_t), target   :: Me
      type(MultigridSolver_t) , pointer :: Child_p          ! Pointer to Child
      integer, intent(in)               :: lvl
      integer, intent(in)               :: nEqn
      !-----------------------------------------------------
      integer                           :: i ! counter

      if (lvl == 1) then
         ! smoothing on the coarsest level
         call Smooth( Me % A, Me % x, Me % b, Me % DimPrb, pre_smooths(lvl) , S_SMOOTHER, Me % BJSmoother)
      else
         Child_p => Me % Child ! set local pointer to a child
         Me % Child % Parent => Me ! set child's parent pointer to this level (Me)

         call Smooth( Me % A, Me % x, Me % b, Me % DimPrb, pre_smooths(lvl) , S_SMOOTHER, Me % BJSmoother) ! Pre smoothing
         ! Compute residual
         Me % r = CSR_MatVecMul( Me % A, Me % x ) 
         do i = 1 , Me % DimPrb
            Me % r(i) = Me % b(i) - Me % r(i)
         end do 

         call MG_Restriction( Me, nEqn, 'solres' ) ! Restrict to a coarse grid
         Me % child % b = Me % child % r
         Me % child % x = 0.0_RP
         call MG_VCycle( Me % Child, lvl-1, nEqn )

         call MG_Prolongation( Me, nEqn )

         call Smooth( Me % A, Me % x, Me % b, Me % DimPrb, pos_smooths(lvl), S_SMOOTHER, Me % BJSmoother) ! Post smoothing
      end if 
   end subroutine MG_VCycle
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine MG_ComputeJacobians(Me,lvl,ComputeTimeDerivative,Time,nEqn)
      use DenseMatUtilities
      implicit none
!-----Arguments---------------------------------------------------
      type(MultigridSolver_t), target    :: Me
      type(MultigridSolver_t) , pointer  :: Child_p          ! Pointer to Child
      procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
      integer, intent(in)                :: nEqn
      integer, intent(in)                :: lvl
!-----Local-Variables---------------------------------------------
      integer             :: i,j
      real(kind=RP)       :: Time
      character(len=1024) :: filename


      call Me % Jacobian % Compute (Me % p_sem, nEqn, time, Me % A, ComputeTimeDerivative)
      select case (S_SMOOTHER)
      case (S_POINTJAC)
         print *, "It's allright"
      case (S_BLOCKJAC)
         call getDBDfromCSR(Me % A,Me % Atmp) ! Get Jacobian diag-blocks
         call ComputeBlockPrec(Me % BJSmoother)


!~         do i = 1, Me % BJSmoother % A_p % num_of_Blocks
!~            ! Visualize Block matrix
!~            write (filename,"(I2.2)") i
!~            filename='mat'//trim(filename)//'.dat'
!~            open(i, file = trim(filename), status = 'new')  
!~            write(i,*) Me % BJSmoother % A_p % Blocks(i) % Matrix   
!~            close(i) 
!~ 
!~            write (filename,"(I2.2)") i
!~            filename='lum'//trim(filename)//'.dat'
!~            open(i+100, file = filename, status = 'new')  
!~            write(i+100,*) Me % BJSmoother % BlockPrec(i) % PLU   
!~            close(i+100) 
!~
!~            write (filename,"(I2.2)") i
!~            filename='lup'//trim(filename)//'.dat'
!~            open(i+200, file = filename, status = 'new')  
!~            write(i+200,*) Me % BJSmoother % BlockPrec(i) % LUpivots   
!~            close(i+200) 
!~         end do
!~
!~         ! Visualization of CSR Jacobian
!~         write (filename,"(I2.2)") lvl
!~         filename='JacFull_lvl'//trim(filename)//'.dat'
!~         call Me % A % Visualize(filename) ! write Jacobian to a file
!~
!~         ! Visualization of Jacobian blocks
!~         call getCSRfromDBD(Me % Atmp,Me % A)
!~         write (filename,"(I2.2)") lvl
!~         filename='JacBlocks_lvl'//trim(filename)//'.dat'
!~         call Me % A % Visualize(filename) ! write Jacobian to a file
      case default
         print *, "Shouldnt be here"
      end select


      ! Compute on a coarser level
      ! *********************************************************
      if (lvl > 1) then
         Child_p => Me % Child ! set local pointer to a child
         Me % Child % Parent => Me ! set child's parent pointer to this level (Me)
         print *, "Calling Jac restriction on lvl: ", lvl
         call MG_Restriction(Me, nEqn, 'soljac')
         ! recursive call
         call MG_ComputeJacobians(Me % Child,lvl-1,ComputeTimeDerivative, Time, nEqn)
      end if
      ! *********************************************************


   end subroutine MG_ComputeJacobians
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_CreateProlongationOperator(this)
      use NodalStorageClass
      implicit none
      !-----------------------------------------------------
      class(MultigridSolver_t), target, intent(inout) :: this
      integer              :: x_Norigin ! Destination polynomial order X
      integer              :: x_Ndest   ! Destination polynomial order X
      type(NodalStorage_t) :: x_spAo    ! Origin nodal storage X
      type(NodalStorage_t) :: x_spAd    ! Destination nodal storage X
      integer              :: y_Norigin ! Destination polynomial order Y
      integer              :: y_Ndest   ! Destination polynomial order Y
      type(NodalStorage_t) :: y_spAo    ! Origin nodal storage Y
      type(NodalStorage_t) :: y_spAd    ! Destination nodal storage Y
      integer              :: z_Norigin ! Destination polynomial order Z
      integer              :: z_Ndest   ! Destination polynomial order Z
      type(NodalStorage_t) :: z_spAo    ! Origin nodal storage Z
      type(NodalStorage_t) :: z_spAd    ! Destination nodal storage Z
      !-----------------------------------------------------

      x_Norigin  = this % Nx
      x_Ndest    = this % parent % Nx
      y_Norigin  = this % Ny
      y_Ndest    = this % parent % Ny
      z_Norigin  = this % Nz
      z_Ndest    = this % parent % Nz
      !-----------------------------------------------------
      call x_spAo % construct(this % p_sem % mesh % nodeType, x_Norigin)
      call x_spAd % construct(this % p_sem % mesh % nodeType, x_Ndest  )
      call y_spAo % construct(this % p_sem % mesh % nodeType, y_Norigin)
      call y_spAd % construct(this % p_sem % mesh % nodeType, y_Ndest  )
      call z_spAo % construct(this % p_sem % mesh % nodeType, z_Norigin)
      call z_spAd % construct(this % p_sem % mesh % nodeType, z_Ndest  )


      ! X directin
      !-----------------------------------------------------
      call PolynomialInterpolationMatrix(x_Norigin, x_Ndest, x_spAo % x, x_spAo % wb, x_spAd % x, this % ProlX)
      !-----------------------------------------------------

      ! Y direction
      !-----------------------------------------------------
      call PolynomialInterpolationMatrix(y_Norigin, y_Ndest, y_spAo % x, y_spAo % wb, y_spAd % x, this % ProlY)
      !-----------------------------------------------------

      ! Z direction
      !-----------------------------------------------------
      call PolynomialInterpolationMatrix(z_Norigin, z_Ndest, z_spAo % x, z_spAo % wb, z_spAd % x, this % ProlZ)
      !-----------------------------------------------------

      ! 3D
      !-----------------------------------------------------
      call MG_Create3DProlongationMatrix(this % Prol3D, x_Norigin, y_Norigin, z_Norigin, x_Ndest, y_Ndest, z_Ndest, &
       x_spAo % x, y_spAo % x, z_spAo % x , x_spAd % x, y_spAd % x, z_spAd % x)
      !-----------------------------------------------------

      call x_spAo % destruct
      call x_spAd % destruct
      call y_spAo % destruct
      call y_spAd % destruct
      call z_spAo % destruct
      call z_spAd % destruct
   end subroutine MG_CreateProlongationOperator
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_CreateRestrictionOperator(this,Rweighted)
      use NodalStorageClass
      implicit none
      !-----------------------------------------------------
      class(MultigridSolver_t), target, intent(inout) :: this
      integer              :: x_Norigin ! x Destination polynomial order
      integer              :: x_Ndest   ! x Destination polynomial order
      type(NodalStorage_t) :: x_spAo    ! x Origin nodal storage
      type(NodalStorage_t) :: x_spAd    ! x Destination nodal storage
      integer              :: y_Norigin ! y Destination polynomial order
      integer              :: y_Ndest   ! y Destination polynomial order
      type(NodalStorage_t) :: y_spAo    ! y Origin nodal storage
      type(NodalStorage_t) :: y_spAd    ! y Destination nodal storage
      integer              :: z_Norigin ! z Destination polynomial order
      integer              :: z_Ndest   ! z Destination polynomial order
      type(NodalStorage_t) :: z_spAo    ! z Origin nodal storage
      type(NodalStorage_t) :: z_spAd    ! z Destination nodal storage
      logical, optional, intent(in) :: Rweighted    ! flag for weighted restriction operator 
      ! only for restriction
      real(kind=RP), allocatable    :: Rtmp(:,:)
      integer                       :: i,j
      !-----------------------------------------------------
      x_Norigin  = this % Nx
      x_Ndest    = this % child % Nx
      y_Norigin  = this % Ny
      y_Ndest    = this % child % Ny
      z_Norigin  = this % Nz
      z_Ndest    = this % child % Nz

      !-----------------------------------------------------
      call x_spAo % construct(this % p_sem % mesh % nodeType, x_Norigin)
      call x_spAd % construct(this % p_sem % mesh % nodeType, x_Ndest  )
      call y_spAo % construct(this % p_sem % mesh % nodeType, y_Norigin)
      call y_spAd % construct(this % p_sem % mesh % nodeType, y_Ndest  )
      call z_spAo % construct(this % p_sem % mesh % nodeType, z_Norigin)
      call z_spAd % construct(this % p_sem % mesh % nodeType, z_Ndest  )

      ! X direction
      !-----------------------------------------------------
      allocate (Rtmp(0:x_Norigin,0:x_Ndest))
      call PolynomialInterpolationMatrix(x_Ndest, x_Norigin, x_spAd % x, x_spAd % wb, x_spAo % x, Rtmp)
      this % RestX = transpose(Rtmp)
      deallocate (Rtmp)
      this % RestJacX = this % RestX 
      if (present(Rweighted) .AND. Rweighted) then
         do j = 0, x_Norigin ; do i = 0, x_Ndest
            this % RestJacX(i,j) = this % RestJacX(i,j) * x_spAo % w(j) / x_spAd % w(i)
         end do            ; end do
      end if
      !-----------------------------------------------------

      ! Y direction
      !-----------------------------------------------------
      allocate (Rtmp(0:y_Norigin,0:y_Ndest))
      call PolynomialInterpolationMatrix(y_Ndest, y_Norigin, y_spAd % x, y_spAd % wb, y_spAo % x, Rtmp)
      this % RestY = transpose(Rtmp)
      deallocate (Rtmp)
      this % RestJacY = this % RestY
      if (present(Rweighted) .AND. Rweighted) then
         do j = 0, y_Norigin ; do i = 0, y_Ndest
            this % RestJacY(i,j) = this % RestJacY(i,j) * y_spAo % w(j) / y_spAd % w(i)
         end do            ; end do
      end if
      !-----------------------------------------------------

      ! Z direction
      !-----------------------------------------------------
      allocate (Rtmp(0:z_Norigin,0:z_Ndest))
      call PolynomialInterpolationMatrix(z_Ndest, z_Norigin, z_spAd % x, z_spAd % wb, z_spAo % x, Rtmp)
      this % RestZ = transpose(Rtmp)
      deallocate (Rtmp)
      this % RestJacZ = this % RestZ
      if (present(Rweighted) .AND. Rweighted) then
         do j = 0, z_Norigin ; do i = 0, z_Ndest
            this % RestJacZ(i,j) = this % RestJacZ(i,j) * z_spAo % w(j) / z_spAd % w(i)
         end do            ; end do
      end if
      !-----------------------------------------------------

      ! 3D
      !-----------------------------------------------------
      call MG_Create3DRestrictionMatrix(this % Rest3D, x_Norigin, y_Norigin, z_Norigin, x_Ndest, y_Ndest, z_Ndest, &
       x_spAo % x, y_spAo % x, z_spAo % x , x_spAd % x, y_spAd % x, z_spAd % x, &
       x_spAo % w, y_spAo % w, z_spAo % w , x_spAd % w, y_spAd % w, z_spAd % w, Rweighted)
      !-----------------------------------------------------

      call x_spAo % destruct
      call x_spAd % destruct
      call y_spAo % destruct
      call y_spAd % destruct
      call z_spAo % destruct
      call z_spAd % destruct

   end subroutine MG_CreateRestrictionOperator
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Prolongation(this,nEqn)
      !
      !---------------------------------------------------------------------
      ! Prolong a 3D array.
      !---------------------------------------------------------------------
      !
      implicit none
      !-----------------------------------------------------
      class(MultigridSolver_t), target, intent(inout) :: this
      class(MultigridSolver_t), pointer               :: child_p
      integer                         , intent(in)    :: nEqn
      integer                                         :: iEl
      integer, dimension(3)                           :: Norigin
      integer, dimension(3)                           :: Ndest
      integer, dimension(nelem)                       :: elsDOF
      !-----------------------------------------------------
      real(kind=RP), dimension(:), allocatable                                          :: x_org
      real(kind=RP), dimension(nEqn,0:this%Nx, 0:this%Ny, 0:this%Nz)                    :: xe_f
      real(kind=RP), dimension(nEqn,0:this%child%Nx, 0:this%child%Ny, 0:this%child%Nz)  :: xe_c
      !-----------------------------------------------------

      child_p => this % child
      Norigin = (/ child_p % Nx, child_p % Ny, child_p % Nz /)
      Ndest = (/ this % Nx, this % Ny, this % Nz /)

      allocate ( x_org(this % DimPrb) ) 
      x_org = this % x 
      xe_f = 0.0_RP

      ! Vec to Elements
      !-----------------------------------------------------
      do iEl = 1, size(child_p % p_sem % mesh % elements)
         elsDOF(iEL) = child_p % p_sem % mesh  % elements(iEL) % storage % NDOF
      end do
      call MG_Vec2El ( child_p % x, child_p % LocalStorage , Norigin , child_p % DimPrb, nEqn, elsDOF , 'x' )
      !-----------------------------------------------------

      do iEl = 1, size(this % p_sem % mesh % elements)
         call MG_Interp3DArrays(nEqn, Norigin, child_p % LocalStorage(iEl) % x , &
                                      Ndest, this % LocalStorage(iEl) % x , &
                                      child_p % ProlX, child_p % ProlY, child_p % ProlZ )
      end do

      ! Elements to Vec
      !-----------------------------------------------------
      do iEl = 1, size(this % p_sem % mesh % elements)
         elsDOF(iEL) = this % p_sem % mesh  % elements(iEL) % storage % NDOF
      end do
      call MG_El2Vec ( this % LocalStorage , this % x, Ndest , this % DimPrb, nEqn, elsDOF , 'x' )
      !-----------------------------------------------------

      this % x = x_org + this % x
      deallocate( x_org )

   end subroutine MG_Prolongation
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Restriction(this, nEqn, r_vars)
      !
      !---------------------------------------------------------------------
      ! Restrict a 3D array..
      !---------------------------------------------------------------------
      !
      implicit none
      !-----------------------------------------------------
      class(MultigridSolver_t), target, intent(inout) :: this
      class(MultigridSolver_t), pointer               :: child_p
      integer                         , intent(in)    :: nEqn
      character(len=*)                , intent(in)    :: r_vars
      integer                                         :: iEl
      integer, dimension(3)                           :: Norigin
      integer, dimension(3)                           :: Ndest
      integer, dimension(nelem)                       :: elsDOF
      integer                                         :: i,j,k ! counters
      !-----------------------------------------------------
      real(kind=RP), dimension(nEqn,0:this%Nx, 0:this%Ny, 0:this%Nz)  :: re_f
      real(kind=RP), dimension(nEqn,0:this%Nx, 0:this%Ny, 0:this%Nz)  :: xe_f
      real(kind=RP), dimension(nEqn,0:this%child%Nx, 0:this%child%Ny, 0:this%child%Nz)  :: re_c
      real(kind=RP), dimension(nEqn,0:this%child%Nx, 0:this%child%Ny, 0:this%child%Nz)  :: xe_c
      !-----------------------------------------------------

      child_p => this % child
      Norigin = (/ this % Nx, this % Ny, this % Nz /)
      Ndest = (/ child_p % Nx, child_p % Ny, child_p % Nz /)

      select case(r_vars)

         case('solres')

            ! Vec to Elements
            !-----------------------------------------------------
            do iEl = 1, size(this % p_sem % mesh % elements)
               elsDOF(iEL) = this % p_sem % mesh  % elements(iEL) % storage % NDOF
            end do
            call MG_Vec2El ( this % x, this % LocalStorage , Norigin , this % DimPrb, nEqn, elsDOF , 'x' )
            call MG_Vec2El ( this % r, this % LocalStorage , Norigin , this % DimPrb, nEqn, elsDOF , 'r' )
            !-----------------------------------------------------

            do iEl = 1, size(this % p_sem % mesh % elements)
               call MG_Interp3DArrays(nEqn, Norigin, this % LocalStorage(iEL) % x , &
                                            Ndest, child_p % LocalStorage(iEl) % x , &
                                            this % RestJacX, this % RestJacY, this % RestJacZ )
                                            !this % RestX, this % RestY, this % RestZ )

               call MG_Interp3DArrays(nEqn, Norigin, this % LocalStorage(iEl) % r , &
                                            Ndest, child_p % LocalStorage(iEl) % r, &
                                            this % RestJacX, this % RestJacY, this % RestJacZ )
                                            !this % RestX, this % RestY, this % RestZ )
            end do

            ! Elements to Vec
            !-----------------------------------------------------
            do iEl = 1, size(child_p % p_sem % mesh % elements)
               elsDOF(iEL) = child_p % p_sem % mesh  % elements(iEL) % storage % NDOF
            end do
            call MG_El2Vec ( child_p % LocalStorage , child_p % x, Ndest , this % DimPrb, nEqn, elsDOF , 'x' )
            call MG_El2Vec ( child_p % LocalStorage , child_p % r, Ndest , this % DimPrb, nEqn, elsDOF , 'r')
            !-----------------------------------------------------

         case('soljac')

            do iEl = 1, size(this % p_sem % mesh % elements)
               call MG_Interp3DArrays(nEqn, Norigin, this % p_sem % mesh % elements(iEl) % storage % Q, &
                                            Ndest, child_p % p_sem % mesh % elements(iEl) % storage % Q, &
                                            this % RestJacX, this % RestJacY, this % RestJacZ )
                                            !this % RestX, this % RestY, this % RestZ )
            end do

      end select

   end subroutine MG_Restriction
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   subroutine MG_Interp3DArrays(Nvars, Nin, inArray, Nout, outArray, InterpX, InterpY, InterpZ )
      implicit none
      !-------------------------------------------------------
      integer                                                        , intent(in)    :: Nvars
      integer      , dimension(3)                                    , intent(in)    :: Nin
      integer      , dimension(3)                                    , intent(in)    :: Nout
      real(kind=RP), dimension(Nvars,0:Nin (1), 0:Nin (2), 0:Nin (3)), intent(in)    :: inArray
      real(kind=RP), dimension(Nvars,0:Nout(1), 0:Nout(2), 0:Nout(3)), intent(out)   :: outArray
      real(kind=RP), dimension(0:Nout(1), 0:Nin(1))                  , intent(in)    :: InterpX
      real(kind=RP), dimension(0:Nout(2), 0:Nin(2))                  , intent(in)    :: InterpY
      real(kind=RP), dimension(0:Nout(3), 0:Nin(3))                  , intent(in)    :: InterpZ
      !-------------------------------------------------------
      integer :: i,j,k,l,m,n
      !-------------------------------------------------------
      
      outArray = 0.0_RP
      
      do n = 0, Nin(3)  ; do k = 0, Nout(3)
         do m = 0, Nin(2)  ; do j = 0, Nout(2)   
            do l = 0, Nin(1)  ; do i = 0, Nout(1)
               outArray(:,i,j,k) = outArray(:,i,j,k) +   InterpX (i,l) &
                                                       * InterpY (j,m) &
                                                       * InterpZ (k,n) &
                                                       * inArray(:,l,m,n)
            end do             ; end do
         end do             ; end do
      end do             ; end do

   end subroutine MG_Interp3DArrays
!
!////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Create3DProlongationMatrix(Mat,N1x,N1y,N1z,N2x,N2y,N2z,x1,y1,z1,x2,y2,z2)
      implicit none
!
!     -----------------------------------------------------------
!     Creates a 3D Lagrange interpolation matrix from a grid with 
!     coordinates x1, y1, z1 (origin) to a grid with coordinates
!     x2, y2, z2 (destination).
!     Original version by : David Kopriva
!     This version by     : Wojciech Laskwoski
!     -----------------------------------------------------------
!
      real(kind=rp)               ,intent(inout) :: Mat((N2x + 1) * (N2y + 1) * (N2z + 1),(N1x + 1) * (N1y + 1) * (N1z + 1)) !<>  
      integer                     ,intent(in)    :: N1x,N1y,N1z  !<  Origin order
      integer                     ,intent(in)    :: N2x,N2y,N2z  !<  Destination order
      real(kind=rp), dimension(:) ,intent(in)    :: x1(0:N1x),y1(0:N1y),z1(0:N1z)     !<  Nodes in origin
      real(kind=rp), dimension(:) ,intent(in)    :: x2(0:N2x),y2(0:N2y),z2(0:N2z)     !<  Nodes in destination
      !----------------------------------------------------------
      integer :: i,j,k,l,m,n      ! Coordinate counters
      integer :: r,s              ! Matrix index counters
      ! integer :: NDOFEL1, NDOFEL2 ! Degrees of freedom in origin and destination
      !----------------------------------------------------------
      
      ! NDOFEL2 = (N2x + 1) * (N2y + 1) * (N2z + 1)
      ! NDOFEL1 = (N1x + 1) * (N1y + 1) * (N1z + 1)
      ! ALLOCATE(Mat(NDOFEL2,NDOFEL1))

      ! print *, "max row: ", N2x + N2y*(N2x + 1) + N2z*(N2x + 1)*(N2y + 1) + 1
      ! print *, "max col: ", N1x + N1y*(N1x + 1) + N1z*(N1x + 1)*(N1y + 1) + 1 
      ! print *, N2x, N2y, N2z
      ! print *, N1x, N1y, N1z

      do k=0, N1z
         do j=0, N1y
            do i=0, N1x
               r = i + j*(N1x + 1) + k*(N1x + 1)*(N1y + 1) + 1            ! Column index
               do n=0, N2z
                  do m=0, N2y
                     do l=0, N2x
                        s = l + m*(N2x + 1) + n*(N2x + 1)*(N2y + 1) + 1   ! Row index
                        
                        ! print *, "i: ", s, " j: ", r, " val: ", MG_LagrangeInterpolationNoBar(x2(l),N1x,x1,i) * &
                        ! MG_LagrangeInterpolationNoBar(y2(m),N1y,y1,j) * MG_LagrangeInterpolationNoBar(z2(n),N1z,z1,k)

                        Mat(s,r) =  MG_LagrangeInterpolationNoBar(x2(l),N1x,x1,i) * &
                                    MG_LagrangeInterpolationNoBar(y2(m),N1y,y1,j) * &
                                    MG_LagrangeInterpolationNoBar(z2(n),N1z,z1,k)
                     end do
                  end do
               end do
            end do
         end do
      end do

   end subroutine MG_Create3DProlongationMatrix
!
!////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Create3DRestrictionMatrix(Mat,N1x,N1y,N1z,N2x,N2y,N2z,x1,y1,z1,x2,y2,z2,w1x,w1y,w1z,w2x,w2y,w2z,Rweighted)
      implicit none
!
!     -----------------------------------------------------------
!     Creates an L2-3D Lagrange interpolation matrix from a grid  
!     with coordinates x1, y1, z1 (origin) to a grid with coordinates
!     x2, y2, z2 (destination)
!     -----------------------------------------------------------
!
      ! real(kind=rp)               ,intent(inout) :: Mat((N1x + 1) * (N1y + 1) * (N1z + 1),(N2x + 1) * (N2y + 1) * (N2z + 1)) !<>
      real(kind=rp)               ,intent(inout) :: Mat((N2x + 1) * (N2y + 1) * (N2z + 1),(N1x + 1) * (N1y + 1) * (N1z + 1)) !<>  
      integer                     ,intent(in)    :: N1x,N1y,N1z  !<  Origin order
      integer                     ,intent(in)    :: N2x,N2y,N2z  !<  Destination order
      real(kind=rp), dimension(:) ,intent(in)    :: x1 (0:N1x),y1 (0:N1y),z1 (0:N1z)     !<  Nodes in origin
      real(kind=rp), dimension(:) ,intent(in)    :: x2 (0:N2x),y2 (0:N2y),z2 (0:N2z)     !<  Nodes in destination
      real(kind=rp), dimension(:) ,intent(in)    :: w1x(0:N1x),w1y(0:N1y),w1z(0:N1z)     !<  Weights in origin
      real(kind=rp), dimension(:) ,intent(in)    :: w2x(0:N2x),w2y(0:N2y),w2z(0:N2z)     !<  Weights in destination
      logical                     ,intent(in)    :: Rweighted !<  flag
      !----------------------------------------------------------
      integer       :: i,j,k,l,m,n      ! Coordinate counters
      integer       :: r,s              ! Matrix index counters
      real(kind=rp) :: MASSterm         ! 
      !----------------------------------------------------------
      
      ! print *, size(Mat,1), size(Mat,2)
      ! print *, Mat(27,125)
      ! error stop "Wojtek"
      ! Create S matrix and store it directly in "Mat"
      do k=0, N1z
         do j=0, N1y
            do i=0, N1x
               r = i + j*(N1x + 1) + k*(N1x + 1)*(N1y + 1) + 1            ! Column index
               do n=0, N2z
                  do m=0, N2y
                     do l=0, N2x
                        s = l + m*(N2x + 1) + n*(N2x + 1)*(N2y + 1) + 1   ! Row index
                        
                        Mat(s,r) = MG_LagrangeInterpolationNoBar(x1(i),N2x,x2,l) * &
                                   MG_LagrangeInterpolationNoBar(y1(j),N2y,y2,m) * &
                                   MG_LagrangeInterpolationNoBar(z1(k),N2z,z2,n) 
                        if ( Rweighted ) Mat(s,r) = Mat(s,r) * w1x(i) * w1y(j) * w1z(k) 

                     end do
                  end do
               end do
            end do
         end do
      end do
      
      if ( Rweighted ) then
         ! Create Mass matrix and finish computing interpolation operator
         do n=0, N2z
            do m=0, N2y
               do l=0, N2x
                  s = l + m*(N2x + 1) + n*(N2x + 1)*(N2y + 1) + 1   ! Row index
      
                  MASSterm = w2x(l) * w2y(m) * w2z(n)

                  ! Matrix Multiplication I = M⁻¹S (taking advantage of the diagonal matrix)
                  Mat(s,:) = Mat(s,:) / MASSterm
               end do
            end do
         end do
      end if
      
   end subroutine MG_Create3DRestrictionMatrix
!
!     ////////////////////////////////////////////////////////////////
!
!     --------------------------------------------------------------------------
!!    Compute the value of the interpolant WITHOUT using the barycentric formula
!     --------------------------------------------------------------------------
!
   FUNCTION MG_LagrangeInterpolationNoBar( x, N, nodes, j) RESULT(l)
!
!     ---------
!     Arguments
!     ---------
!
      use Utilities, only: almostEqual
      REAL(KIND=RP)                 :: l     !>  Lagrange interpolant
      REAL(KIND=RP)                 :: x     !<  Point of evaluation of interpolant
      INTEGER                       :: N     !<  Polynomial order  
      REAL(KIND=RP), DIMENSION(0:N) :: nodes !<  Nodes of Lagrange interpolation
      INTEGER                       :: j     !<  Index of polynomial to be found
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER                       :: i
      REAL(KIND=RP)                 :: numerator, denominator
      REAL(KIND=RP), DIMENSION(0:N) :: values
      !-----------------------------------------------------------------------------
      
      values      = 0.0_RP
      values(j)   = 1.0_RP
      numerator   = 1.0_RP
      denominator = 1.0_RP

      DO i = 0, N
         IF( AlmostEqual( x, nodes(i) ) )    THEN
            l = values(i)
            RETURN 
         ELSE IF (j.ne.i) THEN
         numerator   = numerator*(x - nodes(i))    
         denominator = denominator*(nodes(j) - nodes(i))
         END IF 
      END DO
      l = numerator/denominator

   END FUNCTION MG_LagrangeInterpolationNoBar      
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Vec2El(vec,els,N,vecsize,nEqn,eldof,var2interp)
      implicit none
      !---------------------------------------------------------------------
      integer,       dimension(3),                           intent(in)    :: N
      integer,                                               intent(in)    :: vecsize
      integer,                                               intent(in)    :: nEqn
      integer,       dimension(nelem),                       intent(in)    :: eldof
      real(kind=RP), dimension(vecsize),                     intent(in)    :: vec
      ! real(kind=RP), dimension(nEqn,0:N(1), 0:N(2), 0:N(3)), intent(inout) :: els
      type(TemporaryElementStorage_t) , dimension(nelem)   , intent(inout) :: els
      character(len=*)                                     , intent(in)    :: var2interp
      !---------------------------------------------------------------------
      integer :: eID, firstIdx, lastIdx ! counters
      !---------------------------------------------------------------------
      
      firstIdx = 1
      do eID = 1, nelem
         lastIdx = firstIdx + eldof(eID) * nEqn
         select case (var2interp)
         case ('x')
            els(eID) % x(1:nEqn,0:N(1),0:N(2),0:N(3)) = reshape (vec (firstIdx:lastIdx-1),(/ nEqn, N(1)+1, N(2)+1, N(3)+1/))
         case ('r')
            els(eID) % r(1:nEqn,0:N(1),0:N(2),0:N(3)) = reshape (vec (firstIdx:lastIdx-1),(/ nEqn, N(1)+1, N(2)+1, N(3)+1/))
         end select
         firstIdx = lastIdx
      end do

   end subroutine MG_Vec2El
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_El2Vec(els,vec,N,vecsize,nEqn,eldof,var2interp)
      implicit none
      !---------------------------------------------------------------------
      integer, dimension(3), intent(in)                                    :: N
      integer, intent(in)                                                  :: vecsize
      integer, intent(in)                                                  :: nEqn
      integer, dimension(nelem), intent(in)                                :: eldof
      real(kind=RP), dimension(vecsize), intent(inout)                     :: vec
      ! real(kind=RP), dimension(nEqn,0:N(1), 0:N(2), 0:N(3)), intent(in)     :: els
      type(TemporaryElementStorage_t) , dimension(nelem)   , intent(inout) :: els
      character(len=*)                                     , intent(in)    :: var2interp
      !---------------------------------------------------------------------
      integer           :: eID, firstIdx, lastIdx ! counters
      !---------------------------------------------------------------------
      
      firstIdx = 1
      do eID = 1, nelem
         lastIdx = firstIdx + eldof(eID) * nEqn
         select case (var2interp)
         case ('x')
            vec (firstIdx : lastIdx - 1) = reshape ( els(eID) % x(1:nEqn,0:N(1),0:N(2),0:N(3)) , (/ eldof(eID) * nEqn /) )
         case ('r')
            vec (firstIdx : lastIdx - 1) = reshape ( els(eID) % r(1:nEqn,0:N(1),0:N(2),0:N(3)) , (/ eldof(eID) * nEqn /) )
         end select
         firstIdx = lastIdx
      end do

   end subroutine MG_El2Vec
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
    subroutine MG_SetRHS(this, RHS)
       implicit none
       !-arguments-----------------------------------------------------------
       class(MultigridSolver_t), intent(inout) :: this
       real(kind=RP)           , intent(in)    :: RHS(this % DimPrb)
       !---------------------------------------------------------------------
       
       this % b = RHS
       
    end subroutine MG_SetRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
    function MG_GetX(this) result(x)
       implicit none
       !-----------------------------------------------------------
       class(MultigridSolver_t), intent(inout) :: this
       real(kind=RP)                           :: x(this % DimPrb)
       !-----------------------------------------------------------

       x = this % x

    end function MG_GetX
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function MG_Getxnorm(this,TypeOfNorm) result(xnorm)
      implicit none
      !-----------------------------------------------------------
      class(MultigridSolver_t), intent(inout)  :: this
      character(len=*)                         :: TypeOfNorm
      real(kind=RP)                            :: xnorm
      !-----------------------------------------------------------
      
      select case (TypeOfNorm)
         case ('infinity')
            xnorm = MAXVAL(ABS(this % x))
         case ('l2')
            xnorm = NORM2(this % x)
         case default
            stop 'MultigridSolverClass ERROR: Norm not implemented yet'
      end select

   end function MG_Getxnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function MG_Getrnorm(this) result(rnorm)
      implicit none
!
!     ----------------------------------------
!     Infinity norm
!     ----------------------------------------
!
      !-----------------------------------------------------------
      class(MultigridSolver_t), intent(inout) :: this
      real(kind=RP)                            :: rnorm
      !-----------------------------------------------------------
      real(kind=RP)                            :: residual(this % DimPrb)
      !-----------------------------------------------------------
      
      rnorm = this % rnorm
      
      
   end function MG_Getrnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_SetRHSValue(this, irow, value)
        implicit none
        !-----------------------------------------------------------
        class(MultigridSolver_t) , intent(inout) :: this
        integer                  , intent(in)    :: irow
        real(kind=RP)            , intent(in)    :: value
        !-----------------------------------------------------------
      
        this % b (irow) = value
      
   end subroutine MG_SetRHSValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 !
   subroutine MG_SetRHSValues(this, nvalues, irow, values)
        implicit none
        class(MultigridSolver_t)    , intent(inout)     :: this
        integer                     , intent(in)        :: nvalues
        integer      , dimension(1:), intent(in)        :: irow
        real(kind=RP), dimension(1:), intent(in)        :: values
        !------------------------------------------------------
        integer                                         :: i

        do i=1, nvalues
           if (irow(i)<0) cycle
           this % b(irow(i)) = values(i)
        end do
      
   end subroutine MG_SetRHSValues
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function getnosmooths( inputline, arrdim ) result(n)
      implicit none
!
      integer               :: arrdim, cstart, cend, i, j
      integer, dimension(arrdim) :: n
      character ( len = * ) :: inputline
!
      ! cstart = index(inputline,'[')
      ! cend   = index(inputline, ']', .true. )
      do i = 1, arrdim
         j = 1 + (i-1)*3
         read( inputline( j:j+2 ), * ) n(i)
      end do
!
   end function getnosmooths
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------
!  Convert the DBD matrix to a CSR matrix
!  --------------------------------------
   subroutine getCSRfromDBD(this,Acsr)
      implicit none
      !-arguments---------------------------------------------------
      ! class(Matrix_t), intent(in)    :: this          !<  This matrix
      class(DenseBlockDiagMatrix_t), intent(in)    :: this          !<  This matrix
      class(csrMat_t)              , intent(inout) :: Acsr      !<  Facorized matrix
      !-local-variables---------------------------------------------
      integer :: ii, jj
      integer :: bID
      !-------------------------------------------------------------
      
      if (this % num_of_Rows /= Acsr % num_of_Rows) then
         print*, 'getCSRfromDBD :: ERROR: Matrix dimensions mismatch:', this % num_of_Rows, Acsr % num_of_Rows
         stop
      end if
      
      call Acsr % PreAllocate()
      call Acsr % Reset
      
      call Acsr % SpecifyBlockInfo(this % BlockIdx, this % BlockSizes)
      
!     Fill the Matrix
!     ---------------
      do bID=1, this % num_of_Blocks

         do jj=1, this % BlockSizes(bID)
            do ii=1, this % BlockSizes(bID)
                  call Acsr % SetBlockEntry(bID,bID,ii,jj, this % Blocks (bID) % Matrix(ii,jj))
            end do
         end do

      end do
      
      call Acsr % assembly()
      
   end subroutine getCSRfromDBD
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------
!  Convert the DBD matrix to a CSR matrix
!  Author: Wojciech Laskowski (wj.laskowski@upm.es) 
!  --------------------------------------
   subroutine getDBDfromCSR(this,Adbd)
      implicit none
      !-arguments---------------------------------------------------
      class(csrMat_t),               intent(in)    :: this          
      class(DenseBlockDiagMatrix_t), intent(inout) :: Adbd    
      !-local-variables---------------------------------------------
      integer                     :: i,j,k ! counters
      real(kind=rp), allocatable  :: Mat(:,:) ! local dense matrix (one block/element)
      !-------------------------------------------------------------
      
!     Dimension check
!     ---------------
      if (this % num_of_Rows /= Adbd % num_of_Rows) then
         print *, 'getDBDfromCSR :: ERROR: Matrix dimensions mismatch:', this % num_of_Rows, Adbd % num_of_Rows
         stop
      end if

!     Fill the DBD Matrix
!     ---------------
      if (size(Adbd%BlockSizes,1) .lt. 1) then
         ERROR Stop "BDB Matrix not allocated correctly."
      end if 
      allocate( Mat(Adbd%BlockSizes(1),Adbd%BlockSizes(1))) ! allocate local matrix, this one if for const. pol. order for all elements
      do i = 1, Adbd % num_of_blocks
         Mat = 0.d0
         do j = 1, Adbd % BlockSizes(i) 
            do k = this % Rows(j + (i-1) * Adbd % BlockSizes(i) ), this % Rows(j + (i-1) * Adbd % BlockSizes(i) + 1)-1
               if ( this%Cols(k) .gt. (i) * Adbd % BlockSizes(i)) exit ! check if we exceed the diagonal block
               if ( this%Cols(k) .gt. (i-1) * Adbd % BlockSizes(i) ) then ! check if we are in the right diagonal block
                  Mat( j , this % Cols(k) - (i-1) * Adbd % BlockSizes(i) ) = this % Values (k)
               end if
            end do
         end do
         Adbd % Blocks(i) % Matrix = Mat
      end do
      deallocate(Mat) ! deallocate 
      
   end subroutine getDBDfromCSR
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ComputeBlockPrec(this)
      use DenseMatUtilities
      implicit none
      !-------------------------------------------------------------
      class(BJSmooth_t), target, intent(inout) :: this            !  Iterative solver class
      !-------------------------------------------------------------
      integer :: k      ! Counter
      !-------------------------------------------------------------

      ! print *, "Start test"
      ! print *, "---------------------------"
      ! ! print *, this % A_p % Blocks(1) % Matrix
      ! print *, "---------------------------"
      ! this % A_p % Blocks(1) % Matrix = inverse(this % A_p % Blocks(1) % Matrix )
      ! print *, "---------------------------"
      ! ! print *, this % A_p % Blocks(1) % Matrix
      ! print *, "---------------------------"
      ! print *, "End test"
      ! error stop "TBC"
! !$omp parallel do schedule(runtime)
      do k=1, this % A_p % num_of_Blocks
         ! this % A_p % Blocks(k) % Matrix = inverse( this % A_p % Blocks(k) % Matrix ) ! direct factorziation and save at the same place
         call ComputeLU (A        = this % A_p % Blocks(k) % Matrix, &
                         ALU      = this % BlockPrec(k) % PLU, &
                         LUpivots = this % BlockPrec(k) % LUpivots)
      end do
! !$omp end parallel do
      
   end subroutine ComputeBlockPrec
!   
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
end module MultigridSolverClass

!!! Printing
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!------------------------------------------------------
      ! inquire(file="Q_tot.dat", exist=file_exists)
      ! if (file_exists) then 
      !    open(66, file = 'Q_tot.dat', status = 'old')
      ! else
      !    open(66, file = 'Q_tot.dat', status = 'new')
      ! end if
!------------------------------------------------------
      ! print *, "N1x: ", N1x
      ! print *, "N1y: ", N1y
      ! print *, "N1z: ", N1z

      ! print *, "N2x: ", N2x
      ! print *, "N2y: ", N2y
      ! print *, "N2z: ", N2z

      ! print *, "x1: ", x1
      ! print *, "y1: ", y1
      ! print *, "z1: ", z1

      ! print *, "x2: ", x2
      ! print *, "y2: ", y2
      ! print *, "z2: ", z2
!------------------------------------------------------
!------------------------------------------------------
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////