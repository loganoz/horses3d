!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!      FASMultigridClass.f90
!      Created: 2017-04-XX 10:006:00 +0100 
!      By: Andrés Rueda
!
!      Anisotropic version of the FAS Multigrid Class
!        As is, it is only valid for steady-state cases
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
module AnisFASMultigridClass
   use SMConstants
   use ExplicitMethods
   use DGSEMClass
   use Physics
   use TruncationErrorClass
   use InterpolationMatrices
   use MultigridTypes
   
   implicit none
   
   private
   public AnisFASMultigrid_t
   
   !--------------------------------------------------------------------------------
   ! Type for storing coarsening-direction-wise information for the multigrid scheme
   !--------------------------------------------------------------------------------
   type :: MGStorage_t
      type(DGSem)            , pointer           :: p_sem            ! Pointer to DGSem class variable of current system
      type(DGSem)                                :: tempsem          ! sem used to compute the truncation error using the finest grid solution
      type(MGSolStorage_t), allocatable          :: Var(:)           ! Variables stored element-wise
   end type MGStorage_t
   
   !---------------------------------------------
   ! Main type for AnisFASMultigridClass (THE SOLVER)
   !---------------------------------------------
   type :: AnisFASMultigrid_t
      type(MGStorage_t)                          :: MGStorage(3)          ! Storage of important variables for multigrid. Indices are for coarsening direction (x, y, z)
      type(AnisFASMultigrid_t)   , pointer       :: Child                 ! Next coarser multigrid solver
      type(AnisFASMultigrid_t)   , pointer       :: Parent                ! Next finer multigrid solver
      
      contains
         procedure                                  :: construct
         procedure                                  :: solve
         procedure                                  :: destruct      
   end type AnisFASMultigrid_t
   
!
!  ----------------
!  Module variables
!  ----------------
!
   procedure(SmoothIt_t), pointer :: SmoothIt
   
   !! Parameters
   integer, parameter :: MAX_SWEEPS_DEFAULT = 10000
   
   !! Variables
   integer        :: MGlevels(3)    ! Total number of multigrid levels        
   integer        :: MaxN(3)        ! Maximum polynomial order in every direction
   integer        :: deltaN         !                                         ! TODO: deltaN(3)
   integer        :: nelem          ! Number of elements (this is a p-multigrid implementation)
   integer        :: SweepNumPre    ! Number of sweeps pre-smoothing
   integer        :: SweepNumPost   ! Number of sweeps post-smoothing
   integer        :: SweepNumCoarse ! Number of sweeps on coarsest level
   integer        :: MaxSweeps      ! Maximum number of sweeps in a smoothing process
   logical        :: PostFCycle,PostSmooth ! Post smoothing options
   logical        :: MGOutput       ! Display output?
   logical        :: ManSol         ! Does this case have manufactured solutions?
   logical        :: SmoothFine     ! 
   logical        :: EstimateTE = .FALSE. ! Estimate the truncation error!
   real(kind=RP)  :: SmoothFineFrac ! Fraction that must be smoothed in fine before going to coarser level
   real(kind=RP)  :: cfl            ! Advective cfl number
   real(kind=RP)  :: dcfl           ! Diffusive cfl number
   character(len=LINE_LENGTH)    :: meshFileName
   
!========
 contains
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine construct(this,controlVariables,sem)
      use FTValueDictionaryClass
      implicit none
      !-----------------------------------------------------------
      class(AnisFASMultigrid_t) , intent(inout), TARGET    :: this              !<> Anisotropic FAS multigrid solver to be constructed
      type(FTValueDictionary)  , intent(in), OPTIONAL  :: controlVariables  !<  Input variables
      type(DGSem), TARGET                  , OPTIONAL  :: sem               !<  Fine sem class
      !-----------------------------------------------------------
      integer   :: Dir           ! Direction of coarsening
      integer   :: UserMGlvls    ! User defined number of MG levels
      character(len=LINE_LENGTH)                       :: PostSmoothOptions
      !-----------------------------------------------------------
      
      if (.NOT. PRESENT(sem)) stop 'Fatal error: AnisFASMultigrid needs sem.'
      if (.NOT. PRESENT(controlVariables)) stop 'Fatal error: AnisFASMultigrid needs controlVariables.'
      
!
!     ----------------------------------
!     Read important variables from file
!     ----------------------------------
!
      if (.NOT. controlVariables % containsKey("multigrid levels")) then
         print*, 'Fatal error: "multigrid levels" keyword is needed by the AnisFASMultigrid solver'
         STOP
      end if
      
      UserMGlvls = controlVariables % IntegerValueForKey("multigrid levels")
      
      MaxN(1) = MAXVAL(sem%Nx)
      MaxN(2) = MAXVAL(sem%Ny)
      MaxN(3) = MAXVAL(sem%Nz)
      
      MGlevels(1)  = MIN(MaxN(1),UserMGlvls)
      MGlevels(2)  = MIN(MaxN(2),UserMGlvls)
      MGlevels(3)  = MIN(MaxN(3),UserMGlvls)
      
      if (controlVariables % containsKey("delta n")) then
         deltaN = controlVariables % IntegerValueForKey("delta n")
      else
         deltaN = 1
      end if
      
      if (controlVariables % containsKey("mg sweeps pre" ) .AND. &
          controlVariables % containsKey("mg sweeps post") ) then
         SweepNumPre  = controlVariables % IntegerValueForKey("mg sweeps pre")
         SweepNumPost = controlVariables % IntegerValueForKey("mg sweeps post")
      elseif (controlVariables % containsKey("mg sweeps")) then
         SweepNumPre  = controlVariables % IntegerValueForKey("mg sweeps")
         SweepNumPost = SweepNumPre
      else
         SweepNumPre  = 1
         SweepNumPost = 1
      end if
      
      if (controlVariables % containsKey("mg sweeps coarsest")) then
         SweepNumCoarse = controlVariables % IntegerValueForKey("mg sweeps coarsest")
      else
         SweepNumCoarse = (SweepNumPre + SweepNumPost) / 2
      end if
      
!     Read cfl and dcfl numbers
!     -------------------------
      
      if (controlVariables % containsKey("cfl")) then
         cfl = controlVariables % doublePrecisionValueForKey("cfl")
         if (flowIsNavierStokes) then
            if (controlVariables % containsKey("dcfl")) then
               dcfl       = controlVariables % doublePrecisionValueForKey("dcfl")
            else
               ERROR STOP '"cfl" and "dcfl" keywords must be specified for the FAS integrator'
            end if
         end if
      else
         ERROR STOP '"cfl" keyword must be specified for the FAS integrator'
      end if
      
      select case (controlVariables % StringValueForKey("mg smoother",LINE_LENGTH))
         case('RK3')  ; SmoothIt => TakeRK3Step
         case('SIRK')
            !! SmoothIt => TakeSIRKStep
            error stop ':: SIRK smoother not implemented yet'
         case default 
            write(STD_OUT,*) '"mg smoother" not recognized. Defaulting to RK3.'
            SmoothIt => TakeRK3Step
      end select
      
      PostSmoothOptions = controlVariables % StringValueForKey("postsmooth option",LINE_LENGTH)
      if (trim(PostSmoothOptions) == 'f-cycle') then
         PostFCycle = .true.
      elseif (trim(PostSmoothOptions) == 'smooth') then
         PostSmooth = .true.
      end if
      
      if (controlVariables % containsKey("smooth fine")) then
         SmoothFine = .TRUE.
         SmoothFineFrac = controlVariables % doublePrecisionValueForKey("smooth fine")
      else
         SmoothFine = .FALSE.
      end if
      
      if (controlVariables % containsKey("max mg sweeps")) then
         MaxSweeps = controlVariables % IntegerValueForKey("max mg sweeps")
      else
         MaxSweeps = MAX_SWEEPS_DEFAULT
      end if
      
!
!     -----------------------
!     Update module variables
!     -----------------------
!
      MGOutput       = controlVariables % logicalValueForKey("multigrid output")
      plotInterval   = controlVariables % integerValueForKey("output interval")
      ManSol         = sem % ManufacturedSol
      
      this % MGStorage(1) % p_sem => sem
      this % MGStorage(2) % p_sem => sem
      this % MGStorage(3) % p_sem => sem
      
      nelem = SIZE(sem % mesh % elements)  ! Same for all levels (only p-multigrid)
!
!     -----------------------------
!     Construct linked solvers list
!     -----------------------------
!
      call ConstructLinkedSolvers(this, MAXVAL(MGlevels))
!
!     -----------------------------------------------------------------------
!     Construct interpolation operators and DGSem classes for every subsolver
!     -----------------------------------------------------------------------
!
      meshFileName = controlVariables % stringValueForKey("mesh file name", requestedLength = LINE_LENGTH)
      do Dir = 1, 3
         call ConstructFASInOneDirection(this, MGlevels(Dir), controlVariables, Dir)   ! TODO: change argument to meshFileName
      end do
      
   end subroutine construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------
!  Subroutine for creating the linked list of solvers
!  --------------------------------------------------
   recursive subroutine ConstructLinkedSolvers(Solver, lvl)
      implicit none
      !----------------------------------------------
      type(AnisFASMultigrid_t), TARGET  :: Solver           !<> Current solver
      integer                       :: lvl              !<  Current multigrid level
      !----------------------------------------------
!
!     -------------------------------------
!     Create child if not on coarsest level
!     -------------------------------------
!
      if (lvl > 1) then
         allocate  (Solver % Child)       ! TODO: Remove from here!!!
         Solver % Child % Parent => Solver
         
         call ConstructLinkedSolvers(Solver % Child, lvl - 1)
      end if
      
   end subroutine ConstructLinkedSolvers
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------------------
!  Subroutine for creating interpolation operators, multigrid storage and 
!  DGSem classes for every subsolver in one of the three directions
!  ----------------------------------------------------------------------
   recursive subroutine ConstructFASInOneDirection(Solver, lvl, controlVariables,Dir)
      use FTValueDictionaryClass
      use BoundaryConditionFunctions
      implicit none
      type(AnisFASMultigrid_t), TARGET  :: Solver           !<> Current solver
      integer                       :: lvl              !<  Current multigrid level
      type(FTValueDictionary)       :: controlVariables !<  Control variables (for the construction of coarse sems
      integer                       :: Dir              !<  Coarsening direction
      !----------------------------------------------
      integer, dimension(:), pointer :: N1x,N1y,N1z            !   Order of approximation for every element in current solver
      integer, dimension(:), pointer :: N1
      integer, dimension(nelem,3)    :: N2                     !   Order of approximation for every element in child solver
      integer                        :: i,j,k, iEl             !   Counter
      logical                        :: success                ! Did the creation of sem succeed?
      type(DGSem)          , pointer :: p_sem             ! Pointer to Parent's sem
      type(AnisFASMultigrid_t) , pointer :: Child_p           ! Pointer to Child
      !----------------------------------------------
      !
      integer :: Nxyz(3), fd, l
      !--------------------------
      ! Allocate variable storage
      !--------------------------

      allocate (Solver % MGStorage(Dir) % Var(nelem))
      p_sem => Solver % MGStorage(Dir) % p_sem
      
      ! Define N1x, N1y and N1z according to refinement direction
      N1x => p_sem % Nx
      N1y => p_sem % Ny
      N1z => p_sem % Nz
   
!$omp parallel do
      do k = 1, nelem
         allocate(Solver % MGStorage(Dir) % Var(k) % Q    (0:N1x(k),0:N1y(k),0:N1z(k),N_EQN))
         allocate(Solver % MGStorage(Dir) % Var(k) % E    (0:N1x(k),0:N1y(k),0:N1z(k),N_EQN))
         allocate(Solver % MGStorage(Dir) % Var(k) % S    (0:N1x(k),0:N1y(k),0:N1z(k),N_EQN))
         allocate(Solver % MGStorage(Dir) % Var(k) % Scase(0:N1x(k),0:N1y(k),0:N1z(k),N_EQN))
         
         Solver % MGStorage(Dir) % Var(k) % Scase = 0._RP
      end do   
!$omp end parallel do
!
!     -----------------------------------------------------
!     Fill source term if required (manufactured solutions)
!     -----------------------------------------------------
!
      if (ManSol) then
         do iEl = 1, nelem
            
            do k=0, N1z(iEl)
               do j=0, N1y(iEl)
                  do i=0, N1x(iEl)
                     if (flowIsNavierStokes) then
                        call ManufacturedSolutionSourceNS(p_sem % mesh % elements(iEl) % geom % x(:,i,j,k), &
                                                          0._RP, &
                                                          Solver % MGStorage(Dir) % Var(iEl) % Scase (i,j,k,:)  )
                     else
                        call ManufacturedSolutionSourceEuler(p_sem % mesh % elements(iEl) % geom % x(:,i,j,k), &
                                                             0._RP, &
                                                             Solver % MGStorage(Dir) % Var(iEl) % Scase (i,j,k,:)  )
                     end if
                  end do
               end do
            end do
         end do
      end if
      
      if (lvl > 1) then
         Child_p => Solver % Child
!
!        -----------------------------------------------
!        Allocate restriction and prolongation operators
!        -----------------------------------------------
!
         select case (Dir)
            case(1); N1 => N1x
            case(2); N1 => N1y
            case(3); N1 => N1z
         end select
         
!
!        ---------------------------------------------
!        Create restriction and prolongation operators
!        ---------------------------------------------
!
         ! First we assign the parent's orders to the child
         N2(:,1) = N1x
         N2(:,2) = N1y
         N2(:,3) = N1z
         
         ! Now we create the interpolation operators and change the corresponding orders in the child
         do k=1, nelem
            call CreateInterpolationOperators(N1(k),N2(k,Dir), MaxN(Dir), MGlevels(Dir), lvl-1, DeltaN, Solver % MGStorage(Dir) % p_sem % nodes)
         end do
!
!        -----------------------------------------------------------
!        Create DGSEM class for child in the corresponding direction
!        -----------------------------------------------------------
!
         allocate (Child_p % MGStorage(Dir) % p_sem)
         
         call Child_p % MGStorage(Dir) % p_sem % construct &
                                          (controlVariables  = controlVariables,                                         &
                                           externalState     = Solver % MGStorage(Dir) % p_sem % externalState,          &
                                           externalGradients = Solver % MGStorage(Dir) % p_sem % externalGradients,      &
                                           Nx_ = N2(:,1),    Ny_ = N2(:,2),    Nz_ = N2(:,3),                            &
                                           success = success )
         if (.NOT. success) ERROR STOP "Multigrid: Problem creating coarse solver."
         
         Child_p % MGStorage(Dir) % tempsem = Child_p % MGStorage(Dir) % p_sem
         
         call ConstructFASInOneDirection(Solver % Child, lvl - 1, controlVariables,Dir)
         
      end if
      
   end subroutine ConstructFASInOneDirection
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------
!  Driver of the FAS multigrid solving procedure
!  ---------------------------------------------
   subroutine solve(this,timestep,t,TE)
      implicit none
      class(AnisFASMultigrid_t)        , intent(inout) :: this       !<> The AnisFAS
      integer                          , intent(in)    :: timestep   !<  Current time step
      real(kind=RP)                    , intent(in)    :: t          !<  Current simulation time
      type(TruncationError_t), optional, intent(inout) :: TE(:)      !<> Truncation error for all elements. If present, the multigrid solver also estimates the TE
      !-------------------------------------------------
      character(len=LINE_LENGTH)              :: FileName
      integer                                 :: Dir
      !-------------------------------------------------
      
!
!     -------------------------------------------------
!     Prepare everything for p adaptation (if required)
!     -------------------------------------------------
!
      if (PRESENT(TE)) then
         EstimateTE = .TRUE.
         call InitializeForTauEstimation(TE,this % MGStorage(1) % p_sem,NON_ISOLATED_TE)
      else
         EstimateTE = .FALSE.
      end if
      
      ThisTimeStep = timestep
!
!
!     ---------------------------------------------------------
!     Perform a v-cycle in each direction (x,y,z)
!        (inside, the tau estimation is performed if requested)
!     ---------------------------------------------------------
!
      do Dir = 1, 3
         call FASVCycle(this,t,MGlevels(Dir),Dir,TE)
      end do

!
!     -------------------------------------------
!     Perform p-adaptation and destruct variables
!     -------------------------------------------
!
!~      if (PRESENT(pAdaptator) .AND. pAdaptator % Adapt) then
!~         call pAdaptator % pAdaptTE(this % MGStorage(1) % p_sem,pAdaptator % TE)
!~      end if
      
      !! Finish up
      
      EstimateTE = .FALSE.
   end subroutine solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------
!  Recursive subroutine to perform a v-cycle
!  -----------------------------------------
   recursive subroutine FASVCycle(this,t,lvl,Dir,TE)
      implicit none
      !----------------------------------------------------------------------------
      class(AnisFASMultigrid_t), intent(inout), TARGET :: this    !<  Current level solver
      real(kind=RP)                        :: t       !<  Simulation time
      integer                              :: lvl     !<  Current multigrid level
      integer                              :: Dir     !<  Direction in which multigrid will be performed (x:1, y:2, z:3)
      type(TruncationError_t)              :: TE(:)   !>  Variable containing the truncation error estimation
      !----------------------------------------------------------------------------
      integer                       :: iEl        !Element counter
      type(AnisFASMultigrid_t), pointer :: Child_p        !Pointer to child
      integer                       :: N1(3)          !Polynomial orders in origin solver       (Attention: the origin can be parent or child)
      integer                       :: N2(3)          !Polynomial orders in destination solver  (Attention: the origin can be parent or child)
      real(kind=RP)                 :: dt             !Time variables
      real(kind=RP)                 :: PrevRes
      integer                       :: sweepcount
      type(DGSem)         , pointer :: p_sem          !Pointer to the current sem class
      type(MGSolStorage_t), pointer :: Var(:)         !Pointer to the variable storage class
      type(DGSem)         , pointer :: Childp_sem          !Pointer to the current child's sem class
      type(MGSolStorage_t), pointer :: ChildVar(:)         !Pointer to the child's variable storage class
      integer                       :: NumOfSweeps
      !----------------------------------------------------------------------------
!
!     -----------
!     Definitions
!     -----------
!
      p_sem => this % MGStorage(Dir) % p_sem
      Var   => this % MGStorage(Dir) % Var
!
!     -----------------------
!     Pre-smoothing procedure
!     -----------------------
!
      if (lvl == 1) then
         NumOfSweeps = SweepNumCoarse
      else
         NumOfSweeps = SweepNumPre
      end if
      
      sweepcount = 0
      DO
         do iEl = 1, NumOfSweeps
            dt = MaxTimeStep(p_sem, cfl, dcfl )
            call SmoothIt(p_sem, t, dt )
         end do
         sweepcount = sweepcount + 1
         if (MGOutput) call PlotResiduals( lvl, sweepcount*NumOfSweeps , p_sem )
         
         if (SmoothFine .AND. lvl > 1) then
            call MGRestrictToChild(this,Dir,lvl,t,TE)
            call ComputeTimeDerivative(this % Child % MGStorage(Dir) % p_sem,t)
            
            if (MAXVAL(ComputeMaxResidual(p_sem)) < SmoothFineFrac * MAXVAL(ComputeMaxResidual &
                                                            (this % Child % MGStorage(Dir) % p_sem))) exit
         else
            exit
         end if
      end do
      PrevRes = MAXVAL(ComputeMaxResidual(p_sem))
      
!~       if (MGOutput) call PlotResiduals( lvl , p_sem )
      
      if (lvl > 1) then
         Childp_sem => this % Child % MGStorage(Dir) % p_sem
         ChildVar   => this % Child % MGStorage(Dir) % Var
         
         if (.not. SmoothFine) call MGRestrictToChild(this,Dir,lvl,t,TE)
!
!        --------------------
!        Perform V-Cycle here
!        --------------------
!
         call FASVCycle(this % Child,t, lvl-1, Dir,TE)
!
!        -------------------------------------------
!        Interpolate coarse-grid error to this level
!        -------------------------------------------
!
!$omp parallel
!$omp do private(N1,N2) schedule(runtime)
         do iEl = 1, nelem
            N1 = Childp_sem % mesh % elements(iEl) % Nxyz
            N2 =      p_sem % mesh % elements(iEl) % Nxyz
            
            call Interp3DArraysOneDir(N_EQN, &
                                      N1, ChildVar(iEl) % E, &
                                      N2, Var     (iEl) % E, &
                                      Dir)
         end do
!$omp end do
!$omp barrier
!
!        -----------------------------------------------
!        Correct solution with coarse-grid approximation
!        -----------------------------------------------
!
!$omp do schedule(runtime)
         do iEl = 1, nelem
            p_sem % mesh % elements(iEl) % storage % Q = p_sem % mesh % elements(iEl) % storage % Q + Var(iEl) % E
         end do
!$omp end do
!$omp end parallel
      
      end if
!
!     ------------------------
!     Post-smoothing procedure
!     ------------------------
!
      if (lvl == 1) then
         NumOfSweeps = SweepNumCoarse
      else
         NumOfSweeps = SweepNumPost
      end if
      
      sweepcount=0
      DO
         
         do iEl = 1, NumOfSweeps
            dt = MaxTimeStep(p_sem, cfl, dcfl )
            call SmoothIt(p_sem, t, dt)
         end do

         sweepcount = sweepcount + 1
         if (MGOutput) call PlotResiduals( lvl,sweepcount*NumOfSweeps, p_sem )
         
         if (PostSmooth .or. PostFCycle) then
            if (MAXVAL(ComputeMaxResidual(p_sem)) < PrevRes) exit
         else
            exit
         end if
         
      end do
!~       if (MGOutput) call PlotResiduals( lvl , p_sem )
!
!     -------------------------
!     Compute coarse-grid error
!     -------------------------
!
      if (lvl < MGlevels(Dir)) then
!$omp parallel do schedule(runtime)
         do iEl = 1, nelem
            Var(iEl) % E = p_sem % mesh % elements(iEl) % storage % Q - Var(iEl) % Q
         end do
!$omp end parallel do
      end if
      
   end subroutine FASVCycle
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------
!  Subroutine that restricts to child..... 
!  ------------------------------------------
   subroutine MGRestrictToChild(this,Dir,lvl,t,TE)
      implicit none
      !-------------------------------------------------------------
      class(AnisFASMultigrid_t), TARGET, intent(inout) :: this     !<  Current level solver
      integer                , intent(in)    :: Dir
      integer                , intent(in)    :: lvl
      real(kind=RP)          , intent(in)    :: t
      type(TruncationError_t), intent(inout) :: TE(:)   !>  Variable containing the truncation error estimation 
      !-------------------------------------------------------------
      type(DGSem)          , pointer      :: p_sem          !Pointer to the current sem class
      type(MGSolStorage_t) , pointer      :: Var(:)         !Pointer to the variable storage class
      type(DGSem)          , pointer      :: Childp_sem          !Pointer to the current child's sem class
      type(MGSolStorage_t) , pointer      :: ChildVar(:)         !Pointer to the child's variable storage class
      integer  :: iEl
      integer  :: N1(3)
      integer  :: N2(3)
      !-------------------------------------------------------------
      
      p_sem      => this % MGStorage(Dir) % p_sem
      Var        => this % MGStorage(Dir) % Var
      Childp_sem => this % Child % MGStorage(Dir) % p_sem
      ChildVar   => this % Child % MGStorage(Dir) % Var
      
!$omp parallel
!$omp do private(N1,N2) schedule(runtime)
      do iEl = 1, nelem
         N1 = p_sem      % mesh % elements(iEl) % Nxyz
         N2 = Childp_sem % mesh % elements(iEl) % Nxyz
         
!
!        Restrict solution
!        -----------------
         
         call Interp3DArraysOneDir(N_EQN, &
                                   N1, p_sem      % mesh % elements(iEl) % storage % Q, &
                                   N2, Childp_sem % mesh % elements(iEl) % storage % Q, &
                                   Dir)
!
!        Restrict residual
!        -----------------
            
         call Interp3DArraysOneDir(N_EQN, &
                                   N1, p_sem % mesh % elements(iEl) % storage % Qdot, &
                                   N2, ChildVar(iEl) % S, &
                                   Dir)
      end do
!$omp end do

!
!     **********************************************************************
!     **********************************************************************
!              Now arrange all the storage in the child solver (and estimate TE if necessary...)
!     **********************************************************************
!     **********************************************************************

!
!     ------------------------------------
!     Copy fine grid solution to MGStorage
!        ... and clear source term
!     ------------------------------------
!
!$omp barrier
!$omp do schedule(runtime)
      do iEl = 1, nelem
         ChildVar(iEl) % Q = Childp_sem % mesh % elements(iEl) % storage % Q
         Childp_sem   % mesh % elements(iEl) % storage % S = 0._RP
      end do
!$omp end do
!$omp end parallel
!
!     -------------------------------------------
!     If not on finest level, correct source term
!     -------------------------------------------
!
      if (EstimateTE) call EstimateTruncationError(TE,Childp_sem,t,ChildVar,Dir)
     
     call ComputeTimeDerivative(Childp_sem,t)
      
      
!$omp parallel do schedule(runtime)
      do iEl = 1, nelem
         Childp_sem % mesh % elements(iEl) % storage % S = ChildVar(iEl) % S - Childp_sem % mesh % elements(iEl) % storage % Qdot
      end do
!$omp end parallel do
      
      
   end subroutine MGRestrictToChild
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------
!  Subroutine that destructs an AnisFAS integrator
!  TODO: finish this
!  -----------------------------------------------
   subroutine destruct(this)       
      implicit none
      !-----------------------------------------------------------
      class(AnisFASMultigrid_t), intent(inout) :: this
      !-----------------------------------------------------------
      integer                              :: Dir
      !-----------------------------------------------------------
      
      do Dir=1, 3
         call DestructStorageOneDir(this,MGlevels(Dir),Dir)
      end do
      
      call FinalizeDestruction(this,MAXVAL(MGlevels))
      
   end subroutine destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine DestructStorageOneDir(Solver,lvl,Dir)
      implicit none
      !-----------------------------------------------------------
      class(AnisFASMultigrid_t), intent(inout) :: Solver
      integer                              :: lvl
      integer                              :: Dir
      !-----------------------------------------------------------
      
      !-----------------------------------------------------------
      
      ! First go to finest level (in this direction)
      if (lvl > 1) call DestructStorageOneDir(Solver % Child,lvl-1,Dir)
      
      !Destruct Multigrid storage
      deallocate (Solver % MGStorage(Dir) % Var) ! allocatable components are automatically deallocated
      
      if (lvl < MGlevels(Dir)) then
         call Solver % MGStorage(Dir) % tempsem % destruct()
         call Solver % MGStorage(Dir) % p_sem % destruct()
         deallocate (Solver % MGStorage(Dir) % p_sem)
      else
         nullify    (Solver % MGStorage(Dir) % p_sem)
      end if
      
      
   end subroutine DestructStorageOneDir
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine FinalizeDestruction(Solver,lvl)
      implicit none
      !-----------------------------------------------------------
      class(AnisFASMultigrid_t), intent(inout) :: Solver
      integer                              :: lvl
      !-----------------------------------------------------------
      
      if (lvl > 1) then
         call FinalizeDestruction(Solver % Child,lvl-1)
         deallocate (Solver % Child)
      end if
      
      if (lvl < MAXVAL(MGlevels)) NULLifY(Solver % Parent)
      
   end subroutine FinalizeDestruction
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---
!  Qdot inside an element ignoring the neighbors
!  ---
!~   subroutine IsolatedQdot( self, time, k  )
!~         use DGTimeDerivativeMethods
!~         use ElementClass
!~         implicit none 
!~         !---------------------------------
!~         type(DGSem)   :: self
!~         real(kind=RP) :: time   !<  Time where we compute time derivative
!~         integer       :: k      !<  Element number
!~         !---------------------------------
!~         integer :: Nx, Ny, Nz, NB(2)
!~         integer :: i,j,iFace
!~         real(kind=RP) :: bvExt(N_EQN), flux(N_EQN)
!~         !---------------------------------
!~!
!~!        -----------------------------------------
!~!        Element orders of analyzed element
!~!        -----------------------------------------
!~!
!~         Nx = self%mesh%elements(k)%Nxyz(1)
!~         Ny = self%mesh%elements(k)%Nxyz(2)
!~         Nz = self%mesh%elements(k)%Nxyz(3)
!~!
!~!        -----------------------------------------
!~!        Prolongation of the solution to the faces
!~!        -----------------------------------------
!~!
!~         call ProlongToFaces( self % mesh % elements(k), self % spA(Nx,Ny,Nz) )
!~!
!~!        -------------------------------------------------------
!~!        "Inviscid Riemann fluxes from the solutions on the faces"
!~!        -------------------------------------------------------
!~!
!~         do iFace = 1, 6 
!~            NB = self % mesh % elements(k) % Nxyz (axisMap(:,iFace))
!~            do j = 0, NB(2)
!~               do i = 0, NB(1)
                  
!~                  bvExt = self % mesh % elements(k) % Qb(:,i,j,iFace)
                  
!~                  call RiemannSolver(QLeft  = self % mesh % elements(k) % Qb(:,i,j,iFace), &
!~                                     QRight = bvExt, &
!~                                     nHat   = self % mesh % elements(k) % geom % normal(:,i,j,iFace), &
!~                                     flux   = flux)
!~                  self % mesh % elements(k) % FStarb(:,i,j,iFace) = flux * self % mesh % elements(k) % geom % scal(i,j,iFace)
!~               end do
!~            end do
!~         end do
         
!~         if ( flowIsNavierStokes )     then
!~!
!~!           --------------------------------------
!~!           Set up the face Values on each element
!~!           --------------------------------------
!~!
!~            do iFace = 1, 6 
!~               NB = self % mesh % elements(k) % Nxyz (axisMap(:,iFace))
!~               do j = 0, NB(2)
!~                  do i = 0, NB(1)
!~                     call GradientValuesForQ (Q  = self % mesh % elements(k) % Qb(:,i,j,iFace), &
!~                                              U  = self % mesh % elements(k) % Ub(:,i,j,iFace) )
!~                  end do
!~               end do
!~            end do
!~!
!~!           -----------------------------------
!~!           Compute the gradients over the mesh
!~!           -----------------------------------
!~!
!~            call ComputeDGGradient( self % mesh % elements(k), self % spA(Nx,Ny,Nz), time )
!~!
!~!           ----------------------------------
!~!           Prolong the gradients to the faces
!~!           ----------------------------------
!~!
!~            call ProlongGradientToFaces( self % mesh % elements(k), self % spA(Nx,Ny,Nz) )
!~!
!~!           -------------------------
!~!           Compute gradient averages
!~!           -------------------------
!~!
!~!           This is automatically done by last step

!~         end if

!~!
!~!        ------------------------
!~!        Compute time derivatives
!~!        ------------------------
!~!
!~         call LocalTimeDerivative( self % mesh % elements(k), self % spA(Nx,Ny,Nz), time )

!~   end subroutine IsolatedQdot
   
end module AnisFASMultigridClass