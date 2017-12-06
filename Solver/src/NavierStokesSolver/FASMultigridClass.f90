!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!      FASMultigridClass.f90
!      Created: 2017-04-XX 10:006:00 +0100 
!      By: Andrés Rueda
!
!      FAS Multigrid Class
!        As is, it is only valid for steady-state cases
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
module FASMultigridClass
   use FTValueDictionaryClass
   use SMConstants
   use ExplicitMethods
   use DGSEMClass
   use Physics
   
   use PolynomialInterpAndDerivsModule
   use GaussQuadrature
   implicit none
   
   private
   public FASMultigrid_t
   
   type :: MGStorage_t
      real(kind=RP), dimension(:,:,:,:), allocatable :: Q     ! Solution of the conservative level (before the smoothing)
      real(kind=RP), dimension(:,:,:,:), allocatable :: E     ! Error (for correction)
      real(kind=RP), dimension(:,:,:,:), allocatable :: S     ! Source term interpolated from next finer grid
      real(kind=RP), dimension(:,:,:,:), allocatable :: Scase ! Source term from the specific case that is running (this is actually not necessary for the MG scheme, but it's needed to estimate the truncation error) .. Currently, it only considers the source term from manufactured solutions (state of the code when this module was written)
   end type MGStorage_t
   
   type :: FASMultigrid_t
      type(DGSem)            , pointer           :: p_sem                 ! Pointer to DGSem class variable of current system
      type(FASMultigrid_t)   , pointer           :: Child                 ! Next coarser multigrid solver
      type(FASMultigrid_t)   , pointer           :: Parent                ! Next finer multigrid solver
      integer                                    :: MGlevel               ! Current Multigrid level
      type(MGStorage_t)          , allocatable   :: MGStorage(:)
      
      type(Interpolator_t)       , allocatable   :: Restriction(:,:,:)    ! Restriction operators (element level)
      type(Interpolator_t)       , allocatable   :: Prolongation(:,:,:)   ! Prolongation operators (element level)
      
      contains
         !Subroutines:
         procedure                               :: construct
         procedure                               :: solve
         procedure                               :: destruct
      
   end type FASMultigrid_t
   
   abstract interface
      subroutine SmoothIt_t(sem, t, dt)
         use DGSEMClass
         type(DGSem)   :: sem
         real(kind=RP) :: t, dt
      end subroutine SmoothIt_t
   end interface
   
   procedure(SmoothIt_t), pointer :: SmoothIt
!
!  ----------------
!  Module variables
!  ----------------
!
   !! Parameters !!
   integer, parameter :: MAX_SWEEPS_DEFAULT = 10000
   
   real(kind=RP)  :: cfl
   
   ! Multigrid
   integer        :: MGlevels       ! Total number of multigrid levels
   integer        :: deltaN         ! 
   integer        :: nelem          ! Number of elements
   integer        :: ThisTimeStep   ! Current time step
   integer        :: plotInterval   ! Read to display output
   logical        :: MGOutput       ! Display output?
   logical        :: FMG = .FALSE.  ! Use Full Multigrid algorithm?
   integer        :: SweepNumPre    ! Number of sweeps pre-smoothing
   integer        :: SweepNumPost   ! Number of sweeps post-smoothing
   integer        :: SweepNumCoarse ! Number of sweeps on coarsest level
   integer        :: MaxSweeps      ! Maximum number of sweeps in a smoothing process
   logical        :: PostFCycle,PostSmooth ! Post smoothing options
   !
   logical        :: SmoothFine     !      
   real(kind=RP)  :: SmoothFineFrac ! Fraction that must be smoothed in fine before going to coarser level
   !
   logical        :: ManSol         ! Does this case have manufactured solutions?
   
!========
 contains
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine construct(this,controlVariables,sem)
      implicit none
      !-----------------------------------------------------------
      class(FASMultigrid_t) , intent(inout), target :: this
      type(FTValueDictionary)  , intent(IN), OPTIONAL  :: controlVariables
      type(DGSem), target                  , OPTIONAL  :: sem
      character(len=LINE_LENGTH)                       :: PostSmoothOptions
      !-----------------------------------------------------------
      !Module variables: MGlevels, deltaN
      
      if (.NOT. PRESENT(sem)) stop 'Fatal error: FASMultigrid needs sem.'
      if (.NOT. PRESENT(controlVariables)) stop 'Fatal error: FASMultigrid needs controlVariables.'
      
!
!     ----------------------------------
!     Read important variables from file
!     ----------------------------------
!
      if (.NOT. controlVariables % containsKey("multigrid levels")) then
         print*, 'Fatal error: "multigrid levels" keyword is needed by the FASMultigrid solver'
         STOP
      end if
      
      MGlevels  = controlVariables % IntegerValueForKey("multigrid levels")
      
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
      
      if (controlVariables % containsKey("cfl")) then
         cfl = controlVariables % doublePrecisionValueForKey("cfl")
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
      MGlevels       = MIN (MGlevels,MAX(MAXVAL(sem%Nx),MAXVAL(sem%Ny),MAXVAL(sem%Nz)))
      
      write(STD_OUT,*) 'Constructuing FAS Multigrid'
      write(STD_OUT,*) 'Number of levels:', MGlevels
      
      this % p_sem => sem
      
      nelem = SIZE(sem % mesh % elements)
!
!     --------------------------
!     Create linked solvers list
!     --------------------------
!
      call RecursiveConstructor(this, sem % Nx, sem % Ny, sem % Nz, MGlevels, controlVariables)
      
   end subroutine construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine RecursiveConstructor(Solver, N1x, N1y, N1z, lvl, controlVariables)
      use BoundaryConditionFunctions
      implicit none
      type(FASMultigrid_t), target  :: Solver
      integer, dimension(:)            :: N1x,N1y,N1z      !<  Order of approximation for every element in current solver
      integer                          :: lvl              !<  Current multigrid level
      type(FTValueDictionary)          :: controlVariables !< Control variables (for the construction of coarse sems
      !----------------------------------------------
      integer, dimension(nelem) :: N2x,N2y,N2z            !   Order of approximation for every element in child solver
      integer                   :: N1xMAX,N1yMAX,N1zMAX   !   Maximum polynomial orders for current (fine) grid
      integer                   :: i,j,k, iEl             !   Counter
      logical                   :: success                ! Did the creation of sem succeed?
      type(FASMultigrid_t) , pointer :: Child_p           ! Pointer to Child
      integer                   :: Q1,Q2,Q3,Q4            ! Sizes of vector Q (conserved solution) used for allocation. In this version the argument MOLD of ALLOCATE is not used since several versions of gfortran don't support it yet...
      !----------------------------------------------
      !
      integer :: Nxyz(3), fd, l
      
      Solver % MGlevel = lvl
!
!     --------------------------
!     Allocate Multigrid storage
!     --------------------------
!
      ALLOCATE (Solver % MGStorage(nelem))
!$omp parallel do private(Q1,Q2,Q3,Q4)
      DO k = 1, nelem
         Q1 = SIZE(Solver % p_sem % mesh % elements(k) % storage % Q,1)
         Q2 = SIZE(Solver % p_sem % mesh % elements(k) % storage % Q,2) - 1
         Q3 = SIZE(Solver % p_sem % mesh % elements(k) % storage % Q,3) - 1
         Q4 = SIZE(Solver % p_sem % mesh % elements(k) % storage % Q,4) - 1
         ALLOCATE(Solver % MGStorage(k) % Q    (Q1,0:Q2,0:Q3,0:Q4))
         ALLOCATE(Solver % MGStorage(k) % E    (Q1,0:Q2,0:Q3,0:Q4))
         ALLOCATE(Solver % MGStorage(k) % S    (Q1,0:Q2,0:Q3,0:Q4))
         ALLOCATE(Solver % MGStorage(k) % Scase(Q1,0:Q2,0:Q3,0:Q4))
         
         Solver % MGStorage(k) % Scase = 0._RP
      end DO   
!$omp end parallel do
!
!     --------------------------------------------------------------
!     Fill MGStorage(iEl) % Scase if required (manufactured solutions)
!        (only for lower meshes)
!     --------------------------------------------------------------
!
      if (ManSol) then
         DO iEl = 1, nelem
            
            DO k=0, Solver % p_sem % Nz(iEl)
               DO j=0, Solver % p_sem % Ny(iEl)
                  DO i=0, Solver % p_sem % Nx(iEl)
                     if (flowIsNavierStokes) then
                        call ManufacturedSolutionSourceNS(Solver % p_sem % mesh % elements(iEl) % geom % x(:,i,j,k), &
                                                          0._RP, &
                                                          Solver % MGStorage(iEl) % Scase (:,i,j,k)  )
                     else
                        call ManufacturedSolutionSourceEuler(Solver % p_sem % mesh % elements(iEl) % geom % x(:,i,j,k), &
                                                             0._RP, &
                                                             Solver % MGStorage(iEl) % Scase (:,i,j,k)  )
                     end if
                  end DO
               end DO
            end DO
         end DO
      end if
      
      if (lvl > 1) then
         ALLOCATE  (Solver % Child)
         Child_p => Solver % Child
         Solver % Child % Parent => Solver
!
!        -----------------------------------------------
!        Allocate restriction and prolongation operators
!        -----------------------------------------------
!      
         N1xMAX = MAXVAL(N1x)
         N1yMAX = MAXVAL(N1y)
         N1zMAX = MAXVAL(N1z)
         
         ALLOCATE (Solver  % Restriction (0:N1xMAX,0:N1yMAX,0:N1zMAX))
         ALLOCATE (Child_p % Prolongation(0:N1xMAX,0:N1yMAX,0:N1zMAX))
         
!
!        ---------------------------------------------
!        Create restriction and prolongation operators
!        ---------------------------------------------
!
         DO k=1, nelem
            call CreateInterpolationOperators(Solver % Restriction, Child_p % Prolongation, &
                                              N1x(k),N1y(k),N1z(k),                         &
                                              N2x(k),N2y(k),N2z(k), DeltaN)    ! TODO: Add lobatto flag if required
            
         end DO
         
         ! Create DGSEM class for child
         ALLOCATE (Child_p % p_sem)
         
         call Child_p % p_sem % construct (controlVariables = controlVariables,                                          &
                                           externalState     = Solver % p_sem % externalState,                           &
                                           externalGradients = Solver % p_sem % externalGradients,                       &
                                           Nx_ = N2x,    Ny_ = N2y,    Nz_ = N2z,                                        &
                                           success = success )
         if (.NOT. success) ERROR STOP "Multigrid: Problem creating coarse solver."
         
         call RecursiveConstructor(Solver % Child, N2x, N2y, N2z, lvl - 1, controlVariables)
      end if
      
      
   end subroutine RecursiveConstructor
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------
!  Driver of the FAS multigrid solving procedure
!  ---------------------------------------------
   subroutine solve(this,timestep,t,FullMG,tol)
      implicit none
      class(FASMultigrid_t), intent(inout) :: this
      integer                              :: timestep
      real(kind=RP)                        :: t
      logical           , OPTIONAL         :: FullMG
      real(kind=RP)     , OPTIONAL         :: tol        !<  Tolerance for full multigrid
      !-------------------------------------------------
      integer                                 :: niter
      integer                                 :: i
      character(LEN=LINE_LENGTH)              :: FileName
      !-------------------------------------------------
      
      ThisTimeStep = timestep
      
      if (PRESENT(FullMG) .AND. FullMG) then
         if (.NOT. PRESENT(tol)) ERROR STOP 'FASFMG needs tolerance'
         FMG = .TRUE.
      else
         FMG = .FALSE.
      end if
!
!     -----------------------
!     Perform multigrid cycle
!     -----------------------
!
      if (FMG) then
         call FASFMGCycle(this,t,tol,MGlevels)
      else
         call FASVCycle(this,t,MGlevels,MGlevels)
      end if
      
   end subroutine solve  
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------
!  Recursive subroutine to perform a v-cycle
!  -----------------------------------------
   recursive subroutine FASVCycle(this,t,lvl,MGlevels)
      implicit none
      !----------------------------------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: this     !<  Current level solver
      real(kind=RP)        , intent(in)    :: t        !<  Simulation time
      integer              , intent(in)    :: lvl      !<  Current multigrid level
      integer              , intent(in)    :: MGlevels !<  Number of finest multigrid level
      !----------------------------------------------------------------------------
      integer                       :: iEl,iEQ              !Element/equation counter
      type(FASMultigrid_t), pointer :: Child_p              !Pointer to child
      integer                       :: N1x, N1y, N1z        !Polynomial orders
      integer                       :: N2x, N2y, N2z        !Polynomial orders
      real(kind=RP)                 :: dt                   !Time variables
      real(kind=RP)                 :: maxResidual(N_EQN)
      integer                       :: NumOfSweeps
      real(kind=RP)                 :: PrevRes
      integer                       :: sweepcount           ! Number of sweeps done in a point in time
      !----------------------------------------------------------------------------
      
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
         DO iEl = 1, NumOfSweeps
            dt = MaxTimeStep(this % p_sem, cfl )
            call SmoothIt   (this % p_sem, t, dt )
         end DO
         sweepcount = sweepcount + NumOfSweeps
         
         if (MGOutput) call PlotResiduals( lvl , sweepcount,this % p_sem )
         
         if (SmoothFine .AND. lvl > 1) then ! .AND. .not. FMG
            if (FMG .and. MAXVAL(ComputeMaxResidual(this % p_sem)) < 0.1_RP) exit
            call MGRestrictToChild(this,lvl-1,t)
            call ComputeTimeDerivative(this % Child % p_sem,t)
            
            if (MAXVAL(ComputeMaxResidual(this % p_sem)) < SmoothFineFrac * MAXVAL(ComputeMaxResidual(this % Child % p_sem))) exit
         else
            exit
         end if
         
         if (sweepcount .ge. MaxSweeps) exit
      end DO
      
      PrevRes = MAXVAL(ComputeMaxResidual(this % p_sem))
      
      
      if (lvl > 1) then
         if (.not. SmoothFine) call MGRestrictToChild(this,lvl-1,t)
!
!        --------------------
!        Perform V-Cycle here
!        --------------------
!
         call FASVCycle(this % Child,t, lvl-1,MGlevels)
         
         Child_p => this % Child
!
!        -------------------------------------------
!        Interpolate coarse-grid error to this level
!        -------------------------------------------
!
!$omp parallel do private(iEQ,N1x,N1y,N1z,N2x,N2y,N2z)
         DO iEl = 1, nelem
            DO iEQ = 1, N_EQN
               N1x = Child_p % p_sem % Nx(iEl)
               N1y = Child_p % p_sem % Ny(iEl)
               N1z = Child_p % p_sem % Nz(iEl)
               N2x = this    % p_sem % Nx(iEl)
               N2y = this    % p_sem % Ny(iEl)
               N2z = this    % p_sem % Nz(iEl)
               call Interpolate3D(Q1 = Child_p % MGStorage(iEl) % E(iEQ,:,:,:)      , &
                                  Q2 = this    % MGStorage(iEl) % E(iEQ,:,:,:)      , &
                                  Interp = Child_p % Prolongation(N2x,N2y,N2z) % Mat, &
                                  N1x = N1x,    N1y = N1y,    N1z = N1z             , &
                                  N2x = N2x,    N2y = N2y,    N2z = N2z)
            end DO
         end DO
!$omp end parallel do
!
!        -----------------------------------------------
!        Correct solution with coarse-grid approximation
!        -----------------------------------------------
!
!$omp parallel do
         DO iEl = 1, nelem
            this % p_sem % mesh % elements(iEl) % storage % Q = &
                              this % p_sem % mesh % elements(iEl) % storage % Q + this % MGStorage(iEl) % E
         end DO
!$omp end parallel do
      
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
      
      sweepcount = 0
      DO
         DO iEl = 1, NumOfSweeps
            dt = MaxTimeStep(this % p_sem, cfl )
            call SmoothIt   (this % p_sem, t, dt)
         end DO
         
         sweepcount = sweepcount + NumOfSweeps
         if (MGOutput) call PlotResiduals( lvl, sweepcount , this % p_sem )
         
         if (sweepcount .ge. MaxSweeps) exit
         
         if (lvl > 1 .and. PostFCycle) then
            if (MAXVAL(ComputeMaxResidual(this % p_sem)) > PrevRes) then
               call MGRestrictToChild(this,lvl-1,t)
               call FASVCycle(this,t,lvl-1,lvl)
            else
               exit
            end if
         elseif (PostSmooth .or. PostFCycle) then
            if (MAXVAL(ComputeMaxResidual(this % p_sem)) < PrevRes) exit
         else
            exit
         end if
         
      end DO
      
!
!     -------------------------
!     Compute coarse-grid error
!     -------------------------
!
      if (lvl < MGlevels) then
!$omp parallel do private(iEQ)
         DO iEl = 1, nelem
            this % MGStorage(iEl) % E = this % p_sem % mesh % elements(iEl) % storage % Q - this % MGStorage(iEl) % Q
         end DO
!$omp end parallel do
      end if
      
   end subroutine FASVCycle

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------
!  Recursive subroutine to perform a full multigrid cycle
!  ------------------------------------------------------
   recursive subroutine FASFMGCycle(this,t,tol,lvl)
      implicit none
      !----------------------------------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: this    !<> Current level solver
      real(kind=RP)        , intent(in)    :: t       !<  Simulation time
      real(kind=RP)        , intent(in)    :: tol     !<  Convergence tolerance
      integer              , intent(in)    :: lvl     !<  Current multigrid level
      !----------------------------------------------------------------------------
      integer        :: iEl, iEQ             ! Element and equation counters
      integer        :: N1x, N1y, N1z        ! Origin polynomial orders
      integer        :: N2x, N2y, N2z        ! Destination polynomial orders
      real(kind=RP)  :: dt                   ! Time variables
      real(kind=RP)  :: maxResidual(N_EQN)   ! Maximum residual in each equation
      integer        :: counter              ! Iteration counter
      !----------------------------------------------------------------------------
!
!     ------------------------------------------
!     At the beginning, go to the coarsest level
!        (the initial condition must be passed)
!     ------------------------------------------
!
      if (lvl > 1) then
!$omp parallel do private(iEQ,N1x,N1y,N1z,N2x,N2y,N2z)
         DO iEl = 1, nelem
            DO iEQ = 1, N_EQN
               N1x = this         % p_sem % Nx(iEl)
               N1y = this         % p_sem % Ny(iEl)
               N1z = this         % p_sem % Nz(iEl)
               N2x = this % Child % p_sem % Nx(iEl)
               N2y = this % Child % p_sem % Ny(iEl)
               N2z = this % Child % p_sem % Nz(iEl)
               call Interpolate3D(Q1 = this         % p_sem % mesh % elements(iEl) % storage % Q(iEQ,:,:,:), &
                                  Q2 = this % Child % p_sem % mesh % elements(iEl) % storage % Q(iEQ,:,:,:), &
                                  Interp = this % Restriction(N1x,N1y,N1z) % Mat            , &
                                  N1x = N1x,    N1y = N1y,    N1z = N1z                     , &
                                  N2x = N2x,    N2y = N2y,    N2z = N2z)
            end DO
         end DO
!$omp end parallel do

         call FASFMGCycle(this % Child,t,tol,lvl-1)
      end if
!
!     ----------------------
!     Perform a V-Cycle here
!     ----------------------
!
      counter = 0
      if (lvl > 1 ) then
         DO
            counter = counter + 1
            call FASVCycle(this,t,lvl,lvl)
            maxResidual = ComputeMaxResidual(this % p_sem)
            if (maxval(maxResidual) <= tol) exit
         end DO
      else
         DO
            counter = counter + 1
            dt = MaxTimeStep(this % p_sem, cfl )
            call SmoothIt   (this % p_sem, t, dt )
            maxResidual = ComputeMaxResidual(this % p_sem)
            
            if (MOD(counter,100)==0) call PlotResiduals( lvl ,counter, this % p_sem)
            if (maxval(maxResidual) <= tol) exit
         end DO
      end if
      call PlotResiduals( lvl ,counter, this % p_sem ,.TRUE.)
!
!     --------------------------------------------------
!     If not on finest, Interpolate to next (finer) grid
!     --------------------------------------------------
! 
      if (lvl < MGlevels) then
!$omp parallel do private(iEQ,N1x,N1y,N1z,N2x,N2y,N2z)
         DO iEl = 1, nelem
            DO iEQ = 1, N_EQN
               N1x = this % p_sem % Nx(iEl)
               N1y = this % p_sem % Ny(iEl)
               N1z = this % p_sem % Nz(iEl)
               N2x = this % Parent % p_sem % Nx(iEl)
               N2y = this % Parent % p_sem % Ny(iEl)
               N2z = this % Parent % p_sem % Nz(iEl)
               call Interpolate3D(Q1 = this          % p_sem % mesh % elements(iEl) % storage % Q(iEQ,:,:,:)      , &
                                  Q2 = this % Parent % p_sem % mesh % elements(iEl) % storage % Q(iEQ,:,:,:)      , &
                                  Interp = this % Prolongation(N2x,N2y,N2z) % Mat, &
                                  N1x = N1x,    N1y = N1y,    N1z = N1z             , &
                                  N2x = N2x,    N2y = N2y,    N2z = N2z)
            end DO
         end DO
!$omp end parallel do
      end if
   end subroutine FASFMGCycle

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------
!  Subroutine for restricting the solution and residual to the child solver
!  ------------------------------------------
   subroutine MGRestrictToChild(this,lvl,t)
      implicit none
      !-------------------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: this     !<  Current level solver
      integer              , intent(IN)    :: lvl
      real(kind=RP)        , intent(IN)    :: t
      !-------------------------------------------------------------
      class(FASMultigrid_t), pointer       :: Child_p  ! The child
      integer  :: iEl
      integer  :: iEQ
      integer  :: N1x,N1y,N1z
      integer  :: N2x,N2y,N2z
      !-------------------------------------------------------------
      
      Child_p => this % Child
!
!        -----------------
!        Restrict solution
!        -----------------
!
!$omp parallel do private(iEQ,N1x,N1y,N1z,N2x,N2y,N2z)
         DO iEl = 1, nelem
            DO iEQ = 1, N_EQN
               N1x = this    % p_sem % Nx(iEl)
               N1y = this    % p_sem % Ny(iEl)
               N1z = this    % p_sem % Nz(iEl)
               N2x = Child_p % p_sem % Nx(iEl)
               N2y = Child_p % p_sem % Ny(iEl)
               N2z = Child_p % p_sem % Nz(iEl)
               call Interpolate3D(Q1 = this    % p_sem % mesh % elements(iEl) % storage % Q(iEQ,:,:,:), &
                                  Q2 = Child_p % p_sem % mesh % elements(iEl) % storage % Q(iEQ,:,:,:), &
                                  Interp = this % Restriction(N1x,N1y,N1z) % Mat            , &
                                  N1x = N1x,    N1y = N1y,    N1z = N1z                     , &
                                  N2x = N2x,    N2y = N2y,    N2z = N2z)
            end DO
         end DO
!$omp end parallel do
!
!        -----------------
!        Restrict residual
!        -----------------
!
!$omp parallel do private(iEQ,N1x,N1y,N1z,N2x,N2y,N2z)
         DO iEl = 1, nelem
            DO iEQ = 1, N_EQN
               N1x = this    % p_sem % Nx(iEl)
               N1y = this    % p_sem % Ny(iEl)
               N1z = this    % p_sem % Nz(iEl)
               N2x = Child_p % p_sem % Nx(iEl)
               N2y = Child_p % p_sem % Ny(iEl)
               N2z = Child_p % p_sem % Nz(iEl)
               call Interpolate3D(Q1 = this    % p_sem % mesh % elements(iEl) % storage % Qdot(iEQ,:,:,:), &
                                  Q2 = Child_p % MGStorage(iEl) % S   (iEQ,:,:,:), &
                                  Interp = this % Restriction(N1x,N1y,N1z) % Mat    , &
                                  N1x = N1x,    N1y = N1y,    N1z = N1z             , &
                                  N2x = N2x,    N2y = N2y,    N2z = N2z)
            end DO
         end DO
!$omp end parallel do
!
!     **********************************************************************
!     **********************************************************************
!              Now arrange all the storage in the child solver
!     **********************************************************************
!     **********************************************************************
!

!     ------------------------------------
!     Copy solution from fine grid to MGStorage
!        ... and clear source term
!     ------------------------------------
!
!$omp parallel do private(iEQ)
      DO iEl = 1, nelem
         Child_p % MGStorage(iEl) % Q = Child_p % p_sem % mesh % elements(iEl) % storage % Q
         Child_p % p_sem % mesh % elements(iEl) % storage % S = 0._RP
      end DO
!$omp end parallel do
!
!     -------------------------------------------
!     If not on finest level, correct source term
!     -------------------------------------------
!      
      call ComputeTimeDerivative(Child_p % p_sem,t) 
      
!$omp parallel do
      DO iEl = 1, nelem
         Child_p % p_sem % mesh % elements(iEl) % storage % S = Child_p % MGStorage(iEl) % S - &
                                                                Child_p % p_sem % mesh % elements(iEl) % storage % Qdot
      end DO
!$omp end parallel do
      
      
   end subroutine MGRestrictToChild
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------
!  Subroutine that destructs a FAS integrator
!  ------------------------------------------
   subroutine destruct(this)       
      implicit none
      !-----------------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: this
      !-----------------------------------------------------------
      
      call RecursiveDestructor(this,MGlevels)
      
   end subroutine destruct
   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   recursive subroutine RecursiveDestructor(Solver,lvl)
      implicit none
      !-----------------------------------------------------------
      class(FASMultigrid_t), intent(inout) :: Solver
      integer                              :: lvl
      !-----------------------------------------------------------
      
      ! First go to finest level
      if (lvl > 1) call RecursiveDestructor(Solver % Child,lvl-1)
      
      !Destruct Multigrid storage
      deallocate (Solver % MGStorage) ! allocatable components are automatically deallocated
      
      if (lvl < MGlevels) then
         
         call Solver % p_sem % destruct()
         deallocate (Solver % p_sem)
         
         deallocate (Solver % Prolongation)
         
         nullify    (Solver % Parent)
      else
         nullify    (Solver % p_sem)
      end if
      
      if (lvl > 1) then
         
         deallocate (Solver % Restriction)
         deallocate (Solver % Child)
      end if
      
   end subroutine RecursiveDestructor
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!

  
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  Routines for interpolation procedures (will probably be moved to another module)
!  
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------------------
!  Creates the restriction and prolongation operators for a certain element
!  for multigrid. Takes into account order anisotropy, but the coarse grid 
!  is constructed by reducing the polynomial order uniformly.
!  ------------------------------------------------------------------------
   subroutine CreateInterpolationOperators(Restriction,Prolongation,N1x,N1y,N1z,N2x,N2y,N2z,DeltaN,Lobatto)
      implicit none
      !-----------------------------------------------------
      type(Interpolator_t), target, intent(inout)  :: Restriction (0:,0:,0:)  !>  Restriction operator
      type(Interpolator_t), target, intent(inout)  :: Prolongation(0:,0:,0:)  !>  Prolongation operator
      integer                     , intent(in)     :: N1x, N1y, N1z           !<  Fine grid order(anisotropic) of the element
      integer                     , intent(out)    :: N2x, N2y, N2z           !>  Coarse grid order(anisotropic) of the element
      integer                     , intent(in)     :: DeltaN                  !<  Interval of reduction of polynomial order for coarser level
      logical, OPTIONAL           , intent(in)     :: Lobatto                 !<  Is the quadrature a Legendre-Gauss-Lobatto representation?
      !-----------------------------------------------------
      type(Interpolator_t), pointer :: rest                    ! Pointer to constructed restriction interpolator
      type(Interpolator_t), pointer :: prol                    ! Pointer to constructed prolongation interpolator
      logical                       :: LGL = .FALSE.           ! Is the quadrature a Legendre-Gauss-Lobatto representation? (false is default)
      integer                       :: i,j,k,l,m,n             ! Index counters
      integer                       :: s,r                     ! Row/column counters for operators
      real(kind=RP), allocatable    :: x1 (:), y1 (:), z1 (:)  ! Position of quadrature points on mesh 1
      real(kind=RP), allocatable    :: w1x(:), w1y(:), w1z(:)  ! Weights for quadrature points on mesh 1
      real(kind=RP), allocatable    :: x2 (:), y2 (:), z2 (:)  ! Position of quadrature points on mesh 2
      real(kind=RP), allocatable    :: w2x(:), w2y(:), w2z(:)  ! Weights for quadrature points on mesh 2
      !-----------------------------------------------------
      
      if (PRESENT(Lobatto) .AND. Lobatto) LGL = .TRUE.
!
!     --------------------------------------
!     Compute order of coarse representation
!     --------------------------------------
!
      N2x = N1x - DeltaN
      N2y = N1y - DeltaN
      N2z = N1z - DeltaN
      
      ! The order must be greater or equal to 0 (Legendre-Gauss quadrature) or 1 (Legendre-Gauss-Lobatto)
      if (LGL) then
         if (N2x < 1) N2x = 1
         if (N2y < 1) N2y = 1
         if (N2z < 1) N2z = 1
      else
         if (N2x < 1) N2x = 1       !! The threshold should actually be zero... Using 1 because max eigenvalue subroutine doesn't support 0
         if (N2y < 1) N2y = 1
         if (N2z < 1) N2z = 1
      end if
      
      ! Return if the operators were already created
      if (Restriction(N1x,N1y,N1z) % Created) RETURN
      
      rest => Restriction (N1x,N1y,N1z)
      prol => Prolongation(N1x,N1y,N1z)
!
!     ----------------------------
!     Allocate important variables
!     ----------------------------
!
      !Nodes and weights
      ALLOCATE(x1 (0:N1x), y1 (0:N1y), z1 (0:N1z), &
               w1x(0:N1x), w1y(0:N1y), w1z(0:N1z), &
               x2 (0:N2x), y2 (0:N2y), z2 (0:N2z), &
               w2x(0:N2x), w2y(0:N2y), w2z(0:N2z))
!
!     ------------------------------------------
!     Obtain the quadrature nodes on (1) and (2)
!     ------------------------------------------
!
      if (LGL) then
         call LegendreLobattoNodesAndWeights(N1x, x1, w1x)
         call LegendreLobattoNodesAndWeights(N1y, y1, w1y)
         call LegendreLobattoNodesAndWeights(N1z, z1, w1z)
         call LegendreLobattoNodesAndWeights(N2x, x2, w2x)
         call LegendreLobattoNodesAndWeights(N2y, y2, w2y)
         call LegendreLobattoNodesAndWeights(N2z, z2, w2z)
      else
         call GaussLegendreNodesAndWeights(N1x, x1, w1x)
         call GaussLegendreNodesAndWeights(N1y, y1, w1y)
         call GaussLegendreNodesAndWeights(N1z, z1, w1z)
         call GaussLegendreNodesAndWeights(N2x, x2, w2x)
         call GaussLegendreNodesAndWeights(N2y, y2, w2y)
         call GaussLegendreNodesAndWeights(N2z, z2, w2z)
      end if
!
!     -----------------------------
!     Fill the restriction operator
!     -----------------------------
!
      call Create3DRestrictionMatrix(rest % Mat,N1x,N1y,N1z,N2x,N2y,N2z,x1,y1,z1,x2,y2,z2,w1x,w1y,w1z,w2x,w2y,w2z)
!
!     ------------------------------
!     Fill the prolongation operator
!     ------------------------------
!
      call Create3DInterpolationMatrix(prol % Mat,N2x,N2y,N2z,N1x,N1y,N1z,x2,y2,z2,x1,y1,z1)
      
      ! All done
      rest % Created = .TRUE.
      prol % Created = .TRUE.
      
   end subroutine CreateInterpolationOperators
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------
!  Internal subroutine to print the residuals
!  ------------------------------------------
   subroutine PlotResiduals( lvl, sweeps , sem, white )
      implicit none
      !--------------------------------------------------------
      integer    , intent(in)           :: lvl
      type(DGSem), intent(in)           :: sem
      integer    , intent(in)           :: sweeps
      logical    , intent(in), OPTIONAL :: white
      !--------------------------------------------------------
      real(kind=RP)             :: maxResiduals(N_EQN)
      character(len=5)          :: color1
      character(len=5)          :: color2
      !--------------------------------------------------------
      
      if (PRESENT(white) .AND. white) then
         color1 = achar(27)//'[00m'
      else
         color1 = achar(27)//'[34m'
      end if
      color2 = achar(27)//'[00m'
      
      if( (MOD( ThisTimeStep+1, plotInterval) == 0) .or. (ThisTimeStep .eq. 0) ) then
         maxResiduals = ComputeMaxResidual(sem)
         write(STD_OUT , 110) color1,'FAS lvl', lvl ,"|","it",sweeps,"|", maxResiduals(IRHO) , "|" , maxResiduals(IRHOU) , &
                                 "|", maxResiduals(IRHOV) , "|" , maxResiduals(IRHOW) , "|" , maxResiduals(IRHOE),color2
      end if
      
      110 format (A,A,I3,X,A,X,A,I8,X,A,X,ES10.3,X,A,X,ES10.3,X,A,X,ES10.3,X,A,X,ES10.3,X,A,X,ES10.3,A)
      
   end subroutine PlotResiduals
!
!/////////////////////////////////////////////////////////////////////////////////////////////////

end module FASMultigridClass
