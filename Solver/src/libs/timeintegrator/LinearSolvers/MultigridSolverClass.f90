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
   public MultigridSolver_t, GenericLinSolver_t
   public NUMERICAL_JACOBIAN, ANALYTICAL_JACOBIAN

   type, extends(GenericLinSolver_t) :: MultigridSolver_t
      type(csrMat_t)                             :: A                     ! Matrix to solve
      real(kind=RP), dimension(:), allocatable   :: x                     ! Solution vector
      real(kind=RP), dimension(:), allocatable   :: b                     ! Right hand side
      real(kind=RP), dimension(:), allocatable   :: r                     ! Residual
      real(kind=RP)                              :: rnorm                 ! L2 norm of residual
      real(kind=RP)                              :: tol                   ! Tolerance
      real(kind=RP)                              :: maxiter               ! Max. # iterations
      integer                                    :: MGlevel               ! Current level
      integer                                    :: Nx                    ! Polynomial in X
      integer                                    :: Ny                    ! Polynomial in Y
      integer                                    :: Nz                    ! Polynomial in Z
      type(MultigridSolver_t), pointer           :: Child                 ! Coarser level: MGlevel-1
      type(MultigridSolver_t), pointer           :: Parent                ! Finer level: MGlevel+1
      ! 1D prolongation/restriction operators in each direction 
      real(kind=RP), dimension(:,:), allocatable :: RestX(:,:)             
      real(kind=RP), dimension(:,:), allocatable :: ProlX(:,:)             
      real(kind=RP), dimension(:,:), allocatable :: RestY(:,:)             
      real(kind=RP), dimension(:,:), allocatable :: ProlY(:,:)             
      real(kind=RP), dimension(:,:), allocatable :: RestZ(:,:)             
      real(kind=RP), dimension(:,:), allocatable :: ProlZ(:,:)             
      ! these are only needed for coarse levels
      
   contains
      ! Subroutines:
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
   ! Multigrid
   integer                            :: no_levels    ! Total number of multigrid levels
   integer, dimension(:), allocatable :: MG_levels_x  ! multigrid levels
   integer, dimension(:), allocatable :: MG_levels_y  ! multigrid levels
   integer, dimension(:), allocatable :: MG_levels_z  ! multigrid levels
   integer, dimension(:), allocatable :: pre_smooths  ! pre-smoothing operations
   integer, dimension(:), allocatable :: pos_smooths  ! post-smoothing operations
   integer                            :: nelem        ! Number of elements
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Construct(this,DimPrb,globalDimPrb,nEqn,controlVariables,sem,MatrixShiftFunc)
     implicit none
     !-arguments-----------------------------------------------------------
     class(MultigridSolver_t),  intent(inout), target :: this
     integer                  , intent(in)            :: DimPrb
     integer                  , intent(in)            :: globalDimPrb        
     integer                  , intent(in)            :: nEqn
     type(FTValueDictionary)  , intent(in), optional  :: controlVariables
     type(DGSem), target                  , optional  :: sem
     procedure(MatrixShift_FCN)                       :: MatrixShiftFunc     ! TODO: Make this optional
     procedure(ComputeTimeDerivative_f)               :: ComputeTimeDerivative
     !---------------------------------------------------------------------`
     character(len=LINE_LENGTH)                 :: pc
     integer                                    :: i
     real(kind=RP)                              :: tmp_1

     ! Check dims, MPI and allocate Jacobian
     call this % GenericLinSolver_t % construct(DimPrb, globalDimPrb, nEqn,controlVariables,sem,MatrixShiftFunc)
     MatrixShift => MatrixShiftFunc ! FIXME: ask whether we need this 

!    ***********************************         
!    Set variables from controlVariables
!    ***********************************
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
             read(pc(i:i),'(i)') MG_levels_x(i)
           end do
           ! Levels in Y
           pc = controlVariables % StringValueForKey("define levels y",LINE_LENGTH)
           do i = 1, len_trim(pc)
             read(pc(i:i),'(i)') MG_levels_y(i)
           end do
           ! Levels in Z
           pc = controlVariables % StringValueForKey("define levels z",LINE_LENGTH)
           do i = 1, len_trim(pc)
             read(pc(i:i),'(i)') MG_levels_z(i)
           end do
       else
           ERROR stop ':: FIXME: Need to finish this.'
       end if
       ! *********************************************************

       ! Smoother
       ! *********************************************************
       if ( controlVariables % containsKey("smoother") ) then
         pc = controlVariables % StringValueForKey("smoother",LINE_LENGTH)
       print *, "Smoother: ", pc
       end if

       if ( controlVariables % containsKey("pre smooths") ) then
         pc = controlVariables % StringValueForKey("pre smooths",LINE_LENGTH)
         pre_smooths = getnosmooths(pc,no_levels)
         ! do i = 1, len_trim(pc)
         !   read(pc(i:i),'(i)') pre_smooths(i)
         ! end do
       print *, "No of pre-smoothing operations: ", pre_smooths
       end if

       if ( controlVariables % containsKey("post smooths") ) then
         pc = controlVariables % StringValueForKey("post smooths",LINE_LENGTH)
         pos_smooths = getnosmooths(pc,no_levels)
         ! do i = 1, len_trim(pc)
         !   read(pc(i:i),'(i)') pos_smooths(i)
         ! end do
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
     ! print *, "DOF of the each element: "
     ! do i = 1, size(this % p_sem % mesh % elements) 
     !    print *, i, this % p_sem % mesh % elements(i) % storage % NDOF
     ! end do
     print *, " ************************** "
  
     ! Construct initial variables for the finest level and call recursive constructor
     ! *********************************************************
     this % MGlevel = no_levels
     this % Nx = MG_levels_x(no_levels)
     this % Ny = MG_levels_y(no_levels)
     this % Nz = MG_levels_z(no_levels)
     print *, "Here? 1"
     call this % A % construct(num_of_Rows = DimPrb, withMPI = .false.)
     call this % Jacobian % Configure (sem % mesh, nEqn, this % A)
     print *, "Here? 2"
     ! call this % A % Visualize("JacF.txt") ! Visualize Jacobian
     call MG_Levels_Construct(this,no_levels,controlVariables,nEqn)
     ! *********************************************************

   end subroutine MG_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine MG_Levels_Construct(Me,lvl,controlVariables,nEqn)
     implicit none
     type(MultigridSolver_t), target   :: Me
     type(MultigridSolver_t) , pointer :: Child_p          ! Pointer to Child
     integer, intent(in)               :: nEqn
     integer, intent(in)               :: lvl
     integer                           :: i

     ! character(len=*) :: meshFileName_
     class(FTValueDictionary), intent(in) :: controlVariables                
     !integer, allocatable              :: Nx(:), Ny(:), Nz(:)
     ! integer                           :: polynomialOrder(3)
     logical                           :: success            

     ! for coarse jacobian
     integer                                 :: JacobianComputation = ANALYTICAL_JACOBIAN
     ! integer                                 :: JacobianComputation = NUMERICAL_JACOBIAN

     allocate ( Me % b(Me % DimPrb) ) 
     allocate ( Me % x(Me % DimPrb) ) 
     allocate ( Me % r(Me % DimPrb) ) 
     print *, "Here? 3"

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

       print *, "My interpolation matrices: "
       call MG_CreateRestrictionOperator  ( Me , .true.)
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
     type(MultigridSolver_t) , pointer                :: pMG     
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
      integer                                         :: i, j

      ! print *, "-------------------------------------------"
      ! print *, "-------------------------------------------"
      ! print *, "-------------------------------------------"
      ! print *, "-------------------------------------------"
      ! construct jacobians
      ! pMG => this
      ! do j = no_levels, 2, -1
         ! call pMG % Jacobian % Compute (pMG % p_sem, nEqn, time, pMG % A, ComputeTimeDerivative)
         ! call MG_Restriction(pMG, nEqn, 'soljac')
         ! call pMG % A % Visualize("J1.txt") ! Visualize Jacobian
         ! pMG => this % Child
         ! this % Child % Parent => pMG
         ! print *, "Size: ", pMG % DimPrb
      ! end do
      ! call this % Jacobian % Compute (this % p_sem, nEqn, time, this % A, ComputeTimeDerivative)
      ! call MG_Restriction(this, nEqn, 'soljac')
      ! call this % Child % Jacobian % Compute (this % Child % p_sem, nEqn, time, this % Child % A, ComputeTimeDerivative)
      ! call this % child % A % Visualize("J2.txt") ! Visualize Jacobian
      ! call this % child % child % A % Visualize("J3.txt") ! Visualize Jacobian

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
      ! call this % A % Visualize("Jacobian.txt") ! Visualize Jacobian

      ! open(1, file = 'b.txt', status = 'new')  
      ! write(1,*) this % b  
      ! close(1) 
      this % niter = 0
      this % maxiter = 10
      this % converged = .false.
      this % tol = 1e-12
     
      this % x = 0.0_RP ! setting initial solution
      do while ( (.not. this % converged ) .and. (this % niter .lt. this % maxiter) )
         call MG_VCycle( this, no_levels, nEqn)

         this % r = CSR_MatVecMul( this % A, this % x ) 
         do i = 1 , this % DimPrb
            this % r(i) = this % b(i) - this % r(i)
         end do 
         this % rnorm = norm2(this % r)
         this % niter = this % niter + 1
         ! print *, "Cycle ", this % niter, ". Res: ", this % rnorm
         if (this % rnorm < this % tol) this % converged = .true.
      end do

      ! ERROR stop ':: TBC Wojtek'
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
         call Smooth( Me % A, Me % x, Me % b, Me % DimPrb, pre_smooths(lvl) )
      else
         Child_p => Me % Child ! set local pointer to a child
         Me % Child % Parent => Me ! set child's parent pointer to this level (Me)

         call Smooth( Me % A, Me % x, Me % b, Me % DimPrb, pre_smooths(lvl) ) ! Pre smoothing
         ! Compute residual
         Me % r = CSR_MatVecMul( Me % A, Me % x ) 
         do i = 1 , Me % DimPrb
            Me % r(i) = Me % b(i) - Me % r(i)
         end do 

         call MG_Restriction( Me, nEqn, 'solres' ) ! Restrict to a coarse grid
         Me % child % b = Me % child % r
         call MG_VCycle( Me % Child, lvl-1, nEqn )

         call MG_Prolongation( Me, nEqn )

         call Smooth( Me % A, Me % x, Me % b, Me % DimPrb, pos_smooths(lvl) ) ! Post smoothing
      end if 
   end subroutine MG_VCycle
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine MG_ComputeJacobians(Me,lvl,ComputeTimeDerivative,Time,nEqn)
     implicit none
     type(MultigridSolver_t), target   :: Me
     type(MultigridSolver_t) , pointer :: Child_p          ! Pointer to Child
     procedure(ComputeTimeDerivative_f)              :: ComputeTimeDerivative
     integer, intent(in)               :: nEqn
     integer, intent(in)               :: lvl
     integer                           :: i
     real(kind=RP)                     :: Time

     call Me % Jacobian % Compute (Me % p_sem, nEqn, time, Me % A, ComputeTimeDerivative)
     ! Compute on a coarser level
     ! *********************************************************
     if (lvl > 1) then
      Child_p => Me % Child ! set local pointer to a child
      Me % Child % Parent => Me ! set child's parent pointer to this level (Me)
      call MG_Restriction(Me, nEqn, 'soljac')
      ! recursive call
      call MG_ComputeJacobians(Me % Child,lvl-1,ComputeTimeDerivative, Time, nEqn)
     end if
     ! *********************************************************

     ! Printing
     ! *********************************************************
     print *, " ************************** "
     print *, "We're on level: ", Me%MGlevel
     print *, "Pol. orders: in x/y/z:", Me%Nx, "/", Me%Ny, "/", Me%Nz
     print *, " ************************** "
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
      integer              :: Norigin ! Destination polynomial order
      integer              :: Ndest   ! Destination polynomial order
      type(NodalStorage_t) :: spAo    ! Origin nodal storage
      type(NodalStorage_t) :: spAd    ! Destination nodal storage
      !-----------------------------------------------------

      ! X direction
      !-----------------------------------------------------
      Norigin  = this % Nx
      Ndest    = this % parent % Nx
      call spAo % construct(this % p_sem % mesh % nodeType, Norigin)
      call spAd % construct(this % p_sem % mesh % nodeType, Ndest  )
      call PolynomialInterpolationMatrix(Norigin, Ndest, spAo % x, spAo % wb, spAd % x, this % ProlX)
      call spAo % destruct
      call spAd % destruct
      !-----------------------------------------------------

      ! Y direction
      !-----------------------------------------------------
      Norigin  = this % Ny
      Ndest    = this % parent % Ny
      call spAo % construct(this % p_sem % mesh % nodeType, Norigin)
      call spAd % construct(this % p_sem % mesh % nodeType, Ndest  )
      call PolynomialInterpolationMatrix(Norigin, Ndest, spAo % x, spAo % wb, spAd % x, this % ProlY)
      call spAo % destruct
      call spAd % destruct
      !-----------------------------------------------------

      ! Z direction
      !-----------------------------------------------------
      Norigin  = this % Nz
      Ndest    = this % parent % Nz
      call spAo % construct(this % p_sem % mesh % nodeType, Norigin)
      call spAd % construct(this % p_sem % mesh % nodeType, Ndest  )
      call PolynomialInterpolationMatrix(Norigin, Ndest, spAo % x, spAo % wb, spAd % x, this % ProlZ)
      call spAo % destruct
      call spAd % destruct
      !-----------------------------------------------------
   end subroutine MG_CreateProlongationOperator
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_CreateRestrictionOperator(this,Rweighted)
      use NodalStorageClass
      implicit none
      !-----------------------------------------------------
      class(MultigridSolver_t), target, intent(inout) :: this
      integer              :: Norigin ! Destination polynomial order
      integer              :: Ndest   ! Destination polynomial order
      type(NodalStorage_t) :: spAo    ! Origin nodal storage
      type(NodalStorage_t) :: spAd    ! Destination nodal storage
      logical, optional, intent(in) :: Rweighted    ! flag for weighted restriction operator 
      ! only for restriction
      real(kind=RP), allocatable    :: Rtmp(:,:)
      integer                       :: i,j
      !-----------------------------------------------------

      ! X direction
      !-----------------------------------------------------
      Norigin  = this % Nx
      Ndest    = this % child % Nx
      call spAo % construct(this % p_sem % mesh % nodeType, Norigin)
      call spAd % construct(this % p_sem % mesh % nodeType, Ndest  )
      allocate (Rtmp(0:Norigin,0:Ndest))
      call PolynomialInterpolationMatrix(Ndest, Norigin, spAd % x, spAd % wb, spAo % x, Rtmp)
      this % RestX = transpose(Rtmp)
      deallocate (Rtmp)
      if (present(Rweighted) .AND. Rweighted) then
         do j = 0, Norigin ; do i = 0, Ndest
            this % RestX(i,j) = this % RestX(i,j) * spAo % w(j) / spAd % w(i)
         end do            ; end do
      end if
      call spAo % destruct
      call spAd % destruct
      !-----------------------------------------------------

      ! Y direction
      !-----------------------------------------------------
      Norigin  = this % Ny
      Ndest    = this % child % Ny
      call spAo % construct(this % p_sem % mesh % nodeType, Norigin)
      call spAd % construct(this % p_sem % mesh % nodeType, Ndest  )
      allocate (Rtmp(0:Norigin,0:Ndest))
      call PolynomialInterpolationMatrix(Ndest, Norigin, spAd % x, spAd % wb, spAo % x, Rtmp)
      this % RestY = transpose(Rtmp)
      deallocate (Rtmp)
      if (present(Rweighted) .AND. Rweighted) then
         do j = 0, Norigin ; do i = 0, Ndest
            this % RestY(i,j) = this % RestY(i,j) * spAo % w(j) / spAd % w(i)
         end do            ; end do
      end if
      call spAo % destruct
      call spAd % destruct
      !-----------------------------------------------------

      ! Z direction
      !-----------------------------------------------------
      Norigin  = this % Nz
      Ndest    = this % child % Nz
      call spAo % construct(this % p_sem % mesh % nodeType, Norigin)
      call spAd % construct(this % p_sem % mesh % nodeType, Ndest  )
      allocate (Rtmp(0:Norigin,0:Ndest))
      call PolynomialInterpolationMatrix(Ndest, Norigin, spAd % x, spAd % wb, spAo % x, Rtmp)
      this % RestZ = transpose(Rtmp)
      deallocate (Rtmp)
      if (present(Rweighted) .AND. Rweighted) then
         do j = 0, Norigin ; do i = 0, Ndest
            this % RestZ(i,j) = this % RestZ(i,j) * spAo % w(j) / spAd % w(i)
         end do            ; end do
      end if
      call spAo % destruct
      call spAd % destruct
      !-----------------------------------------------------

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
      call MG_Vec2El ( child_p % x, xe_c, Norigin , child_p % DimPrb, nEqn, elsDOF )
      !-----------------------------------------------------

      do iEl = 1, size(this % p_sem % mesh % elements)
         call MG_Interp3DArrays(nEqn, Norigin, xe_c, &
                                      Ndest, xe_f, &
                                      child_p % ProlX, child_p % ProlY, child_p % ProlZ )
      end do

      ! Elements to Vec
      !-----------------------------------------------------
      do iEl = 1, size(this % p_sem % mesh % elements)
         elsDOF(iEL) = this % p_sem % mesh  % elements(iEL) % storage % NDOF
      end do
      call MG_El2Vec ( xe_f, this % x, Ndest , this % DimPrb, nEqn, elsDOF )
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
            call MG_Vec2El ( this % x, xe_f, Norigin , this % DimPrb, nEqn, elsDOF )
            call MG_Vec2El ( this % r, re_f, Norigin , this % DimPrb, nEqn, elsDOF )
            !-----------------------------------------------------

            do iEl = 1, size(this % p_sem % mesh % elements)
               call MG_Interp3DArrays(nEqn, Norigin, xe_f, &
                                            Ndest, xe_c, &
                                            this % RestX, this % RestY, this % RestZ )

               call MG_Interp3DArrays(nEqn, Norigin, re_f, &
                                            Ndest, re_c, &
                                            this % RestX, this % RestY, this % RestZ )
            end do

            ! Elements to Vec
            !-----------------------------------------------------
            do iEl = 1, size(child_p % p_sem % mesh % elements)
               elsDOF(iEL) = child_p % p_sem % mesh  % elements(iEL) % storage % NDOF
            end do
            call MG_El2Vec ( xe_c, child_p % x, Ndest , this % DimPrb, nEqn, elsDOF )
            call MG_El2Vec ( re_c, child_p % r, Ndest , this % DimPrb, nEqn, elsDOF )
            !-----------------------------------------------------

         case('soljac')

            do iEl = 1, size(this % p_sem % mesh % elements)
               call MG_Interp3DArrays(nEqn, Norigin, this % p_sem % mesh % elements(iEl) % storage % Q, &
                                            Ndest, child_p % p_sem % mesh % elements(iEl) % storage % Q, &
                                            this % RestX, this % RestY, this % RestZ )
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
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_Vec2El(vec,els,N,vecsize,nEqn,eldof)
      implicit none
      !---------------------------------------------------------------------
      integer,       dimension(3),                           intent(in)    :: N
      integer,                                               intent(in)    :: vecsize
      integer,                                               intent(in)    :: nEqn
      integer,       dimension(nelem),                       intent(in)    :: eldof
      real(kind=RP), dimension(vecsize),                     intent(in)    :: vec
      real(kind=RP), dimension(nEqn,0:N(1), 0:N(2), 0:N(3)), intent(inout) :: els
      !---------------------------------------------------------------------
      integer :: eID, firstIdx, lastIdx ! counters
      !---------------------------------------------------------------------
      
      firstIdx = 1
      do eID = 1, nelem
         lastIdx = firstIdx + eldof(eID) * nEqn
         els(1:nEqn,0:N(1),0:N(2),0:N(3)) = reshape (vec (firstIdx:lastIdx-1),(/ nEqn, N(1)+1, N(2)+1, N(3)+1/))
         firstIdx = lastIdx
      end do

   end subroutine MG_Vec2El
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MG_El2Vec(els,vec,N,vecsize,nEqn,eldof)
      implicit none
      !---------------------------------------------------------------------
      integer, dimension(3), intent(in)                                     :: N
      integer, intent(in)                                                   :: vecsize
      integer, intent(in)                                                   :: nEqn
      integer, dimension(nelem), intent(in)                                 :: eldof
      real(kind=RP), dimension(vecsize), intent(inout)                      :: vec
      real(kind=RP), dimension(nEqn,0:N(1), 0:N(2), 0:N(3)), intent(in)     :: els
      !---------------------------------------------------------------------
      integer           :: eID, firstIdx, lastIdx ! counters
      !---------------------------------------------------------------------
      
      firstIdx = 1
      do eID = 1, nelem
         lastIdx = firstIdx + eldof(eID) * nEqn
         vec (firstIdx : lastIdx - 1) = reshape ( els , (/ eldof(eID) * nEqn /) )
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
end module MultigridSolverClass

!!! OLD STUFF
      !!! class(Matrix_t), allocatable               :: A                     ! Jacobian matrix
      !!! ! print *, sem % mesh % storage % Q 
      !!! allocate(SparseBlockDiagMatrix_t :: this % A)
      !!! print *, "No of elements: ", nelem
      !!! print *, "Ok 1"
      !!! call this % A % construct (num_of_Blocks = nelem)
      !!! print *, "Ok 2"
      !!! call this % Jacobian % Configure (sem % mesh, nEqn, this % A)
      !!! !call this % Jacobian % Compute (sem, nEqn, time, this % A, ComputeTimeDerivative, BlockDiagonalized = .TRUE.)
      !!! print *, "No of blocks: ", this % A % num_of_Blocks
      !!! select type(Mat => this % A)
      !!!   type is(SparseBlockDiagMatrix_t)
      !!!   do i = 1, Mat % num_of_Blocks
      !!!     print *, Mat % Blocks(i) % Matrix % Cols(:)
      !!!   end do
      !!! end select

!!! Printing
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!------------------------------------------------------
! print *, "Rest X: ", Me  % RestX
! print *, "Rest X size: ", size(Me  % RestX,1), size(Me  % RestX,2)
! print *, "Rest Y: ", Me  % RestY
! print *, "Rest Y size: ", size(Me  % RestY,1), size(Me  % RestY,2)
! print *, "Rest Z: ", Me  % RestZ
! print *, "Rest Z size: ", size(Me  % RestZ,1), size(Me  % RestZ,2)
! print *, "Prol X: ", Child_p  % ProlX
! print *, "Prol X size: ", size(Child_p  % ProlX,1), size(Child_p  % ProlX,2)
! print *, "Prol Y: ", Child_p  % ProlY
! print *, "Prol Y size: ", size(Child_p  % ProlY,1), size(Child_p  % ProlY,2)
! print *, "Prol Z: ", Child_p  % ProlZ
! print *, "Prol Z size: ", size(Child_p  % ProlZ,1), size(Child_p  % ProlZ,2)
!------------------------------------------------------
!! just rho:
!print *, "-----------------------------------------------------------"
!print *, " ----- RHO -----"
!print *, "Fine"                                      
!print *, this % p_sem % mesh % elements(iEl) % storage % Q(1,:,:,:)
!print *, "Coarse"                                      
!print *, child_p % p_sem % mesh % elements(iEl) % storage % Q(1,:,:,:)
!print *, " ----- RHO U -----"
!print *, "Fine"                                      
!print *, this % p_sem % mesh % elements(iEl) % storage % Q(2,:,:,:)
!print *, "Coarse"                                      
!print *, child_p % p_sem % mesh % elements(iEl) % storage % Q(2,:,:,:)
!print *, " ----- RHO V -----"
!print *, "Fine"                                      
!print *, this % p_sem % mesh % elements(iEl) % storage % Q(3,:,:,:)
!print *, "Coarse"                                      
!print *, child_p % p_sem % mesh % elements(iEl) % storage % Q(3,:,:,:)
!print *, " ----- RHO W -----"
!print *, "Fine"                                      
!print *, this % p_sem % mesh % elements(iEl) % storage % Q(4,:,:,:)
!print *, "Coarse"                                      
!print *, child_p % p_sem % mesh % elements(iEl) % storage % Q(4,:,:,:)
!print *, " ----- E -----"
!print *, "Fine"                                      
!print *, this % p_sem % mesh % elements(iEl) % storage % Q(5,:,:,:)
!print *, "Coarse"                                      
!print *, child_p % p_sem % mesh % elements(iEl) % storage % Q(5,:,:,:)
!print *, "-----------------------------------------------------------"
!------------------------------------------------------
!       print *, "N_x fine: ", Nin(1), "N_x coarse: ", Nout(1)
!       print *, "N_y fine: ", Nin(2), "N_y coarse: ", Nout(2)
!       print *, "N_z fine: ", Nin(3), "N_z coarse: ", Nout(3)
! 
!       do l = 0, Nin(1)  ; do i = 0, Nout(1)
!          print *, "x :: ", i
!          print *, "y :: ", l
!          print *, "R(x,y) :: ", InterpX (i,l)
!       end do             ; end do
!------------------------------------------------------
!       print *, "We're here 8."
!       print *, " R_x :: ", this % RestX, "size: ", size(this % RestX,1), " x ", size(this % RestX,2)
!       print *, " R_y :: ", this % RestZ, "size: ", size(this % RestY,1), " x ", size(this % RestY,2)
!       print *, " R_z :: ", this % RestY, "size: ", size(this % RestZ,1), " x ", size(this % RestZ,2)
!------------------------------------------------------
      ! print *, "Interesting: ", size(Me % A % Rows,1)
!------------------------------------------------------
! call MG_Interp3DArrays(nEqn, Norigin, this % p_sem % mesh % elements(iEl) % storage % QDot, &
!                              Ndest, child_p % p_sem % mesh % elements(iEl) % storage % QDot, &
!                              this % RestX, this % RestY, this % RestZ )
! call MG_Interp3DArrays(nEqn, Norigin, this % p_sem % mesh % elements(iEl) % storage % U_x, &
!                              Ndest, child_p % p_sem % mesh % elements(iEl) % storage % U_x, &
!                              this % RestX, this % RestY, this % RestZ )
! call MG_Interp3DArrays(nEqn, Norigin, this % p_sem % mesh % elements(iEl) % storage % U_y, &
!                              Ndest, child_p % p_sem % mesh % elements(iEl) % storage % U_y, &
!                              this % RestX, this % RestY, this % RestZ )
! call MG_Interp3DArrays(nEqn, Norigin, this % p_sem % mesh % elements(iEl) % storage % U_z, &
!                              Ndest, child_p % p_sem % mesh % elements(iEl) % storage % U_z, &
!                              this % RestX, this % RestY, this % RestZ )
!------------------------------------------------------
         !print *, "We're here loop over els:", iEl
         !print *, "Org Q: ", iEl                                      
         !print *, this % p_sem % mesh % elements(iEl) % storage % Q
         ! call MG_Interp3DArrays(nEqn, Norigin, this % p_sem % mesh % elements(iEl) % storage % Q, &
         !                              Ndest, child_p % p_sem % mesh % elements(iEl) % storage % Q, &
         !                              this % RestX, this % RestY, this % RestZ )
         ! print *, "Fine x", iEl                                      
         ! print *, "rho: "                                      
         ! print *, xe_f(1,:,:,:)
         ! print *, "rhou: "                                      
         ! print *, xe_f(2,:,:,:)
         ! print *, "rhov: "                                      
         ! print *, xe_f(3,:,:,:)
         ! print *, "rhow: "                                      
         ! print *, xe_f(4,:,:,:)
         ! print *, "E: "                                      
         ! print *, xe_f(5,:,:,:)
         ! print *, "Coarse x after", iEl                                      
         ! print *, "rho: "                                      
         ! print *, xe_c(1,:,:,:)
         ! print *, "rhou: "                                      
         ! print *, xe_c(2,:,:,:)
         ! print *, "rhov: "                                      
         ! print *, xe_c(3,:,:,:)
         ! print *, "rhow: "                                      
         ! print *, xe_c(4,:,:,:)
         ! print *, "E: "                                      
         ! print *, xe_c(5,:,:,:)
!------------------------------------------------------
      ! print *, "Sizes: "
      ! print *, size(this % x)
      ! print *, size(this % r)
      ! print *, size(this % b)
      ! print *, size(this % child % x)
      ! print *, size(this % child % r)
      ! print *, size(this % child % b)
!------------------------------------------------------
!------------------------------------------------------
      ! print *, "Values (Jac. Fine)"
      ! print *, "------------------------------------------------------------------------------"
      ! print *, this % A % Values(:)
      ! print *, "------------------------------------------------------------------------------"
      ! call this % A % Visualize("JacF.txt") ! write Jacobian to a file
      !print *, "Cols"
      !print *, "------------------------------------------------------------------------------"
      !print *, this % A % Cols(:)
      !print *, "------------------------------------------------------------------------------"
      !print *, "Rows"
      !print *, "------------------------------------------------------------------------------"
      !print *, this % A % Rows(:)
      !print *, "------------------------------------------------------------------------------"
      !tmpsize = sizeof(this % A % Values(:)) + sizeof(this % A % Cols(:)) + sizeof(this % A % Rows(:))
      !tmpsize = sizeof(this % A % Values(:)) / sizeof(tol)
      ! print *, "Sizeof Jacobian: ", tmpsize
      ! print *, "Sizeof sem: ", sizeof(this % p_sem % mesh)
!------------------------------------------------------
      ! print *, "-------------------------------------------"
      ! print *, "xe_c : ", xe_c
      ! print *, "-------------------------------------------"
      ! print *, "xe_f : ", xe_f
      ! print *, "-------------------------------------------"
      ! print *, "P_x : ", child_p % ProlX
      ! print *, "P_y : ", child_p % ProlY
      ! print *, "P_z : ", child_p % ProlZ
      ! print *, "-------------------------------------------"
!------------------------------------------------------
      ! print *, "Fine Q:"
      ! print *, "------------------------------------------------------------------------------"
      ! print *, this % p_sem % mesh % elements(1) % storage % Q(1,:,:,1)  
      ! print *, "------------------------------------------------------------------------------"

      ! print *, "Coarse Q:"
      ! print *, "------------------------------------------------------------------------------"
      ! print *, this % child % p_sem % mesh % elements(1) % storage % Q(1,:,:,1)  
      ! print *, "------------------------------------------------------------------------------"

      ! print *, "Coarse Q after Restriction:"
      ! print *, "------------------------------------------------------------------------------"
      ! print *, this % child % p_sem % mesh % storage % Q  
      ! print *, "------------------------------------------------------------------------------"

      ! call this % Child % A % Visualize("JacC_pre.txt") ! write Jacobian to a file

      ! print *, "Values (Jac. Coarse)"
      ! print *, "------------------------------------------------------------------------------"
      ! print *, this % Child % A % Values(:)
      ! print *, "------------------------------------------------------------------------------"
      ! call this % Child % A % Visualize("JacC.txt") ! write Jacobian to a file
!------------------------------------------------------
     ! print *, "x before: ", this % x
     ! print *, "x after: ", this % x
     ! print *, "RHS: ", this % b
     ! print *, "x: ", this % x
!------------------------------------------------------
     ! do j = 1, 5
     !    call Smooth( this % A, this % x, this % b, this % DimPrb, pre_smooths(no_levels) )
     !    this % r = CSR_MatVecMul( this % A, this % x ) ! CSR matrix product
     !    do i=1,this%DimPrb
     !       this % r(i) = this % b(i) - this % r(i)
     !    end do 

     !    this % rnorm = norm2( this % r )
     !    print *, "Residual norm is: ", this % rnorm

     !    call MG_Restriction(this, nEqn, 'solres')

     !    this % child % b = this % child % r
     !    call Smooth( this % child % A, this % child % x, this % child % b, this % child % DimPrb, pre_smooths(no_levels-1) )

     !    call MG_Prolongation(this, nEqn)

     !    call Smooth( this % A, this % x, this % b, this % DimPrb, pos_smooths(no_levels) )
     ! end do
!------------------------------------------------------
!------------------------------------------------------
!------------------------------------------------------
!------------------------------------------------------
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////