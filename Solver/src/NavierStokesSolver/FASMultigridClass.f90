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
   IMPLICIT NONE
   
   PRIVATE
   PUBLIC FASMultigrid_t
   
   TYPE :: MGStorage_t
      REAL(KIND=RP), DIMENSION(:,:,:,:), ALLOCATABLE :: Q     ! Solution of the conservative level (before the smoothing)
      REAL(KIND=RP), DIMENSION(:,:,:,:), ALLOCATABLE :: E     ! Error (for correction)
      REAL(KIND=RP), DIMENSION(:,:,:,:), ALLOCATABLE :: S     ! Source term interpolated from next finer grid
      REAL(KIND=RP), DIMENSION(:,:,:,:), ALLOCATABLE :: Scase ! Source term from the specific case that is running (this is actually not necessary for the MG scheme, but it's needed to estimate the truncation error) .. Currently, it only considers the source term from manufactured solutions (state of the code when this module was written)
   END TYPE MGStorage_t
   
   TYPE :: FASMultigrid_t
      TYPE(DGSem)            , POINTER           :: p_sem                              ! Pointer to DGSem class variable of current system
      TYPE(DGSem)                                :: tempsem   ! sem used to compute the truncation error using the finest grid solution
      ! Variables that are specially needed for Multigrid
      TYPE(FASMultigrid_t)   , POINTER           :: Child                 ! Next coarser multigrid solver
      TYPE(FASMultigrid_t)   , POINTER           :: Parent                ! Next finer multigrid solver
      INTEGER                                    :: MGlevel               ! Current Multigrid level
      TYPE(MGStorage_t)          , ALLOCATABLE   :: MGStorage(:)
      
      TYPE(Interpolator_t)       , ALLOCATABLE   :: Restriction(:,:,:)    ! Restriction operators (element level)
      TYPE(Interpolator_t)       , ALLOCATABLE   :: Prolongation(:,:,:)   ! Prolongation operators (element level)
      
   CONTAINS
      !Subroutines:
      PROCEDURE                                  :: construct
      PROCEDURE                                  :: solve
      PROCEDURE                                  :: destruct
      
   END TYPE FASMultigrid_t
   
!
!  ----------------
!  Module variables
!  ----------------
!
   REAL(KIND=RP)  :: cfl
   CHARACTER(len=LINE_LENGTH) :: SmootherName
   
   ! Multigrid
   INTEGER        :: MGlevels       ! Total number of multigrid levels
   INTEGER        :: deltaN         ! 
   INTEGER        :: nelem          ! Number of elements (this is a p-multigrid implementation)
   INTEGER        :: ThisTimeStep   ! Current time step
   INTEGER        :: plotInterval   ! Read to display output
   LOGICAL        :: MGOutput       ! Display output?
   LOGICAL        :: FMG = .FALSE.  ! Use Full Multigrid algorithm?
   INTEGER        :: SweepNumPre    ! Number of sweeps pre-smoothing
   INTEGER        :: SweepNumPost   ! Number of sweeps post-smoothing
   INTEGER        :: SweepNumC      ! Number of sweeps on coarsest level
   LOGICAL        :: PostFCycle,PostSmooth ! Post smoothing options
   REAL(KIND=RP)  :: SmoothFine     ! Fraction that must be smoothed in fine before going to coarser level
   ! 
   LOGICAL        :: ManSol         ! Does this case have manufactured solutions?
   
!========
 CONTAINS
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE construct(this,controlVariables,sem)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(FASMultigrid_t) , INTENT(INOUT), TARGET :: this
      TYPE(FTValueDictionary)  , INTENT(IN), OPTIONAL  :: controlVariables
      TYPE(DGSem), TARGET                  , OPTIONAL  :: sem
      CHARACTER(len=LINE_LENGTH)                       :: PostSmoothOptions
      !-----------------------------------------------------------
      !Module variables: MGlevels, deltaN
      
      IF (.NOT. PRESENT(sem)) stop 'Fatal error: FASMultigrid needs sem.'
      IF (.NOT. PRESENT(controlVariables)) stop 'Fatal error: FASMultigrid needs controlVariables.'
      
!
!     ----------------------------------
!     Read important variables from file
!     ----------------------------------
!
      IF (.NOT. controlVariables % containsKey("multigrid levels")) THEN
         print*, 'Fatal error: "multigrid levels" keyword is needed by the FASMultigrid solver'
         STOP
      END IF
      
      MGlevels  = controlVariables % IntegerValueForKey("multigrid levels")
      
      IF (controlVariables % containsKey("delta n")) THEN
         deltaN = controlVariables % IntegerValueForKey("delta n")
      ELSE
         deltaN = 1
      END IF
      
      IF (controlVariables % containsKey("mg sweeps pre" ) .AND. &
          controlVariables % containsKey("mg sweeps post") ) THEN
         SweepNumPre  = controlVariables % IntegerValueForKey("mg sweeps pre")
         SweepNumPost = controlVariables % IntegerValueForKey("mg sweeps post")
      ELSEIF (controlVariables % containsKey("mg sweeps")) THEN
         SweepNumPre  = controlVariables % IntegerValueForKey("mg sweeps")
         SweepNumPost = controlVariables % IntegerValueForKey("mg sweeps")
      ELSE
         SweepNumPre  = 1
         SweepNumPost = 1
      END IF
      
      IF (controlVariables % containsKey("mg sweeps coarsest")) THEN
         SweepNumC = controlVariables % IntegerValueForKey("mg sweeps coarsest")
      ELSE
         SweepNumC = (SweepNumPre + SweepNumPost) / 2
      END IF
      
      IF (controlVariables % containsKey("cfl")) THEN
         cfl = controlVariables % doublePrecisionValueForKey("cfl")
      ELSE
         ERROR STOP '"cfl" keyword must be specified for the FAS integrator'
      END IF
      
      SmootherName  = controlVariables % StringValueForKey("mg smoother",LINE_LENGTH)
      
      PostSmoothOptions = controlVariables % StringValueForKey("postsmooth option",LINE_LENGTH)
      if (trim(PostSmoothOptions) == 'f-cycle') then
         PostFCycle = .true.
      elseif (trim(PostSmoothOptions) == 'smooth') then
         PostSmooth = .true.
      end if
      
      IF (controlVariables % containsKey("smooth fine")) THEN
         SmoothFine = controlVariables % doublePrecisionValueForKey("smooth fine")
      ELSE
         SmoothFine = -1._RP
      END IF
      
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
      CALL RecursiveConstructor(this, sem % Nx, sem % Ny, sem % Nz, MGlevels, controlVariables)
      
   END SUBROUTINE construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   RECURSIVE SUBROUTINE RecursiveConstructor(Solver, N1x, N1y, N1z, lvl, controlVariables)
      use BoundaryConditionFunctions
      IMPLICIT NONE
      TYPE(FASMultigrid_t), TARGET  :: Solver
      INTEGER, DIMENSION(:)            :: N1x,N1y,N1z      !<  Order of approximation for every element in current solver
      INTEGER                          :: lvl              !<  Current multigrid level
      TYPE(FTValueDictionary)          :: controlVariables !< Control variables (for the construction of coarse sems
      !----------------------------------------------
      INTEGER, DIMENSION(nelem) :: N2x,N2y,N2z            !   Order of approximation for every element in child solver
      INTEGER                   :: N1xMAX,N1yMAX,N1zMAX   !   Maximum polynomial orders for current (fine) grid
      INTEGER                   :: i,j,k, iEl             !   Counter
      LOGICAL                   :: success                ! Did the creation of sem succeed?
      TYPE(FASMultigrid_t) , POINTER :: Child_p           ! Pointer to Child
      INTEGER                   :: Q1,Q2,Q3,Q4            ! Sizes of vector Q (conserved solution) used for allocation. In this version the argument MOLD of ALLOCATE is not used since several versions of gfortran don't support it yet...
      !----------------------------------------------
      !
      INTEGER :: Nxyz(3), fd, l
      
      Solver % MGlevel = lvl
!
!     --------------------------
!     Allocate Multigrid storage
!     --------------------------
!
      ALLOCATE (Solver % MGStorage(nelem))
!$omp parallel do private(Q1,Q2,Q3,Q4)
      DO k = 1, nelem
         Q1 = SIZE(Solver % p_sem % mesh % elements(k) % storage % Q,1) - 1 
         Q2 = SIZE(Solver % p_sem % mesh % elements(k) % storage % Q,2) - 1
         Q3 = SIZE(Solver % p_sem % mesh % elements(k) % storage % Q,3) - 1
         Q4 = SIZE(Solver % p_sem % mesh % elements(k) % storage % Q,4)
         ALLOCATE(Solver % MGStorage(k) % Q    (0:Q1,0:Q2,0:Q3,Q4))
         ALLOCATE(Solver % MGStorage(k) % E    (0:Q1,0:Q2,0:Q3,Q4))
         ALLOCATE(Solver % MGStorage(k) % S    (0:Q1,0:Q2,0:Q3,Q4))
         ALLOCATE(Solver % MGStorage(k) % Scase(0:Q1,0:Q2,0:Q3,Q4))
         
         Solver % MGStorage(k) % Scase = 0._RP
      END DO   
!$omp end parallel do
!
!     --------------------------------------------------------------
!     Fill MGStorage(iEl) % Scase if required (manufactured solutions)
!        (only for lower meshes)
!     --------------------------------------------------------------
!
      IF (ManSol) THEN
         DO iEl = 1, nelem
            
            DO k=0, Solver % p_sem % Nz(iEl)
               DO j=0, Solver % p_sem % Ny(iEl)
                  DO i=0, Solver % p_sem % Nx(iEl)
                     IF (flowIsNavierStokes) THEN
                        CALL ManufacturedSolutionSourceNS(Solver % p_sem % mesh % elements(iEl) % geom % x(:,i,j,k), &
                                                          0._RP, &
                                                          Solver % MGStorage(iEl) % Scase (i,j,k,:)  )
                     ELSE
                        CALL ManufacturedSolutionSourceEuler(Solver % p_sem % mesh % elements(iEl) % geom % x(:,i,j,k), &
                                                             0._RP, &
                                                             Solver % MGStorage(iEl) % Scase (i,j,k,:)  )
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END IF
      
      IF (lvl > 1) THEN
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
            CALL CreateInterpolationOperators(Solver % Restriction, Child_p % Prolongation, &
                                              N1x(k),N1y(k),N1z(k),                         &
                                              N2x(k),N2y(k),N2z(k), DeltaN)    ! TODO: Add lobatto flag if required
            
         END DO
         
         ! Create DGSEM class for child
         ALLOCATE (Child_p % p_sem)
         
         CALL Child_p % p_sem % construct (controlVariables = controlVariables,                                          &
                                           externalState     = Solver % p_sem % externalState,                           &
                                           externalGradients = Solver % p_sem % externalGradients,                       &
                                           Nx_ = N2x,    Ny_ = N2y,    Nz_ = N2z,                                        &
                                           success = success )
         IF (.NOT. success) ERROR STOP "Multigrid: Problem creating coarse solver."
         
         Child_p % tempsem = Child_p % p_sem
         
         CALL RecursiveConstructor(Solver % Child, N2x, N2y, N2z, lvl - 1, controlVariables)
      END IF
      
      
   END SUBROUTINE RecursiveConstructor
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------
!  Driver of the FAS multigrid solving procedure
!  ---------------------------------------------
   SUBROUTINE solve(this,timestep,t,FullMG,tol)
      IMPLICIT NONE
      CLASS(FASMultigrid_t), INTENT(INOUT) :: this
      INTEGER                              :: timestep
      REAL(KIND=RP)                        :: t
      LOGICAL           , OPTIONAL         :: FullMG
      REAL(KIND=RP)     , OPTIONAL         :: tol        !<  Tolerance for full multigrid
      !-------------------------------------------------
      INTEGER                                 :: niter
      INTEGER                                 :: i
      CHARACTER(LEN=LINE_LENGTH)              :: FileName
      !-------------------------------------------------
      
      ThisTimeStep = timestep
      
      IF (PRESENT(FullMG) .AND. FullMG) THEN
         IF (.NOT. PRESENT(tol)) ERROR STOP 'FASFMG needs tolerance'
         FMG = .TRUE.
      ELSE
         FMG = .FALSE.
      END IF
!
!     -----------------------
!     Perform multigrid cycle
!     -----------------------
!
      IF (FMG) THEN
         CALL FASFMGCycle(this,t,tol,MGlevels)
      ELSE
         CALL FASVCycle(this,t,MGlevels,MGlevels)
      END IF
      
   END SUBROUTINE solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   RECURSIVE SUBROUTINE InterpolateFinest(this,lvl)
      IMPLICIT NONE
      !---------------------------------------------
      CLASS(FASMultigrid_t), INTENT(INOUT) :: this
      INTEGER                              :: lvl
      !---------------------------------------------
      TYPE(FASMultigrid_t), POINTER :: Child_p        !Pointer to child
      INTEGER :: iEl, iEQ, N1x,N1y,N1z, N2x,N2y,N2z
      !---------------------------------------------
      INTEGER :: i,j,k
      
      IF (lvl == MGlevels) this % tempsem = this % p_sem
      
      IF (lvl>1) THEN
!
!        -----------------
!        Restrict solution
!        -----------------
!
         Child_p => this % Child
!$omp parallel do private(iEQ,N1x,N1y,N1z,N2x,N2y,N2z)
         DO iEl = 1, nelem
            DO iEQ = 1, N_EQN
               N1x = this    % p_sem % Nx(iEl)
               N1y = this    % p_sem % Ny(iEl)
               N1z = this    % p_sem % Nz(iEl)
               N2x = Child_p % p_sem % Nx(iEl)
               N2y = Child_p % p_sem % Ny(iEl)
               N2z = Child_p % p_sem % Nz(iEl)
               CALL Interpolate3D(Q1 = this    % tempsem % mesh % elements(iEl) % storage % Q(:,:,:,iEQ), &
                                  Q2 = Child_p % tempsem % mesh % elements(iEl) % storage % Q(:,:,:,iEQ), &
                                  Interp = this % Restriction(N1x,N1y,N1z) % Mat            , &
                                  N1x = N1x,    N1y = N1y,    N1z = N1z                     , &
                                  N2x = N2x,    N2y = N2y,    N2z = N2z)
               
            END DO
            
!~             DO k=0, Child_p % p_sem % Nz(iEl)
!~                DO j=0, Child_p % p_sem % Ny(iEl)
!~                   DO i=0, Child_p % p_sem % Nx(iEl)
!~                      CALL ManufacturedSolutionState(  Child_p % tempsem % mesh % elements(iEl) % geom % x(:,i,j,k), &
!~                                                       0._RP, &
!~                                                       Child_p % tempsem % mesh % elements(iEl) % Q(i,j,k,:))
!~                   END DO
!~                END DO
!~             END DO
!~             IF (iEl == 1) THEN
!~                print*, '------------------------------------------------------------'
!~                print*, N1x, N1y, N1z, N2x, N2y, N2z
!~                print*, Child_p % tempsem % mesh % elements(iEl) % Q(:,:,:,1)
!~                read(*,*)
!~             END IF
         END DO
!$omp end parallel do
      
         CALL InterpolateFinest(this % Child,lvl-1)
      END IF
      
   END SUBROUTINE InterpolateFinest   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------
!  Recursive subroutine to perform a v-cycle
!  -----------------------------------------
   RECURSIVE SUBROUTINE FASVCycle(this,t,lvl,MGlevels)
      IMPLICIT NONE
      !----------------------------------------------------------------------------
      CLASS(FASMultigrid_t), INTENT(INOUT) :: this     !<  Current level solver
      REAL(KIND=RP)        , intent(in)    :: t        !<  Simulation time
      INTEGER              , intent(in)    :: lvl      !<  Current multigrid level
      INTEGER              , intent(in)    :: MGlevels !<  Number of finest multigrid level
      !----------------------------------------------------------------------------
      INTEGER                       :: iEl,iEQ              !Element/equation counter
      TYPE(FASMultigrid_t), POINTER :: Child_p              !Pointer to child
      INTEGER                       :: N1x, N1y, N1z        !Polynomial orders
      INTEGER                       :: N2x, N2y, N2z        !Polynomial orders
      REAL(KIND=RP)                 :: dt                   !Time variables
      REAL(KIND=RP)                 :: maxResidual(N_EQN)
      INTEGER                       :: NumOfSweeps
      real(kind=RP)                 :: PrevRes
      integer                       :: sweepcount
      !----------------------------------------------------------------------------
      
!
!     -----------------------
!     Pre-smoothing procedure
!     -----------------------
!
      IF (lvl == 1) THEN
         NumOfSweeps = SweepNumC
      ELSE
         NumOfSweeps = SweepNumPre
      END IF
      
      sweepcount = 0
      DO
         DO iEl = 1, NumOfSweeps
            dt = MaxTimeStep(this % p_sem, cfl )
            CALL SmoothIt   (this % p_sem, t, dt )
         END DO
         sweepcount = sweepcount + 1
         
         IF (MGOutput) CALL PlotResiduals( lvl , NumOfSweeps*sweepcount,this % p_sem )
         
         IF (SmoothFine > 0._RP .AND. lvl > 1) THEN ! .AND. .not. FMG
            IF (FMG .and. MAXVAL(ComputeMaxResidual(this % p_sem)) < 0.1_RP) EXIT
            CALL MGRestrictToChild(this,lvl-1,t)
            CALL ComputeTimeDerivative(this % Child % p_sem,t)
            
            IF (MAXVAL(ComputeMaxResidual(this % p_sem)) < SmoothFine * MAXVAL(ComputeMaxResidual(this % Child % p_sem))) EXIT
         ELSE
            EXIT
         END IF
      END DO
      
      PrevRes = MAXVAL(ComputeMaxResidual(this % p_sem))
      
      
      IF (lvl > 1) THEN
         IF (SmoothFine < 0._RP) CALL MGRestrictToChild(this,lvl-1,t)
!
!        --------------------
!        Perform V-Cycle here
!        --------------------
!
         CALL FASVCycle(this % Child,t, lvl-1,MGlevels)
         
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
               CALL Interpolate3D(Q1 = Child_p % MGStorage(iEl) % E(:,:,:,iEQ)      , &
                                  Q2 = this    % MGStorage(iEl) % E(:,:,:,iEQ)      , &
                                  Interp = Child_p % Prolongation(N2x,N2y,N2z) % Mat, &
                                  N1x = N1x,    N1y = N1y,    N1z = N1z             , &
                                  N2x = N2x,    N2y = N2y,    N2z = N2z)
            END DO
         END DO
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
         END DO
!$omp end parallel do
      
      END IF
!
!     ------------------------
!     Post-smoothing procedure
!     ------------------------
!
      IF (lvl == 1) THEN
         NumOfSweeps = SweepNumC
      ELSE
         NumOfSweeps = SweepNumPost
      END IF
      sweepcount = 0
      DO
         DO iEl = 1, NumOfSweeps
            dt = MaxTimeStep(this % p_sem, cfl )
            CALL SmoothIt   (this % p_sem, t, dt)
         END DO
         
         sweepcount = sweepcount + 1
         IF (MGOutput) CALL PlotResiduals( lvl, NumOfSweeps*sweepcount , this % p_sem )
         
         IF (lvl > 1 .and. PostFCycle) THEN
            IF (MAXVAL(ComputeMaxResidual(this % p_sem)) > PrevRes) THEN
               CALL MGRestrictToChild(this,lvl-1,t)
               CALL FASVCycle(this,t,lvl-1,lvl)
            ELSE
               EXIT
            END IF
         ELSEIF (PostSmooth .or. PostFCycle) THEN
            IF (MAXVAL(ComputeMaxResidual(this % p_sem)) < PrevRes) exit
         ELSE
            EXIT
         END IF
         
      END DO
      
!
!     -------------------------
!     Compute coarse-grid error
!     -------------------------
!
      IF (lvl < MGlevels) THEN
!$omp parallel do private(iEQ)
         DO iEl = 1, nelem
            this % MGStorage(iEl) % E = this % p_sem % mesh % elements(iEl) % storage % Q - this % MGStorage(iEl) % Q
         END DO
!$omp end parallel do
      END IF
      
   END SUBROUTINE FASVCycle

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------
!  Recursive subroutine to perform a full multigrid cycle
!  ------------------------------------------------------
   RECURSIVE SUBROUTINE FASFMGCycle(this,t,tol,lvl)
      IMPLICIT NONE
      !----------------------------------------------------------------------------
      CLASS(FASMultigrid_t), INTENT(INOUT) :: this    !<> Current level solver
      REAL(KIND=RP)        , intent(in)    :: t       !<  Simulation time
      REAL(KIND=RP)        , intent(in)    :: tol     !<  Convergence tolerance
      INTEGER              , intent(in)    :: lvl     !<  Current multigrid level
      !----------------------------------------------------------------------------
      INTEGER        :: iEl, iEQ             ! Element and equation counters
      INTEGER        :: N1x, N1y, N1z        ! Origin polynomial orders
      INTEGER        :: N2x, N2y, N2z        ! Destination polynomial orders
      REAL(KIND=RP)  :: dt                   ! Time variables
      REAL(KIND=RP)  :: maxResidual(N_EQN)   ! Maximum residual in each equation
      INTEGER        :: counter              ! Iteration counter
      !----------------------------------------------------------------------------
!
!     ------------------------------------------
!     At the beginning, go to the coarsest level
!        (the initial condition must be passed)
!     ------------------------------------------
!
      IF (lvl > 1) THEN
!$omp parallel do private(iEQ,N1x,N1y,N1z,N2x,N2y,N2z)
         DO iEl = 1, nelem
            DO iEQ = 1, N_EQN
               N1x = this         % p_sem % Nx(iEl)
               N1y = this         % p_sem % Ny(iEl)
               N1z = this         % p_sem % Nz(iEl)
               N2x = this % Child % p_sem % Nx(iEl)
               N2y = this % Child % p_sem % Ny(iEl)
               N2z = this % Child % p_sem % Nz(iEl)
               CALL Interpolate3D(Q1 = this         % p_sem % mesh % elements(iEl) % storage % Q(:,:,:,iEQ), &
                                  Q2 = this % Child % p_sem % mesh % elements(iEl) % storage % Q(:,:,:,iEQ), &
                                  Interp = this % Restriction(N1x,N1y,N1z) % Mat            , &
                                  N1x = N1x,    N1y = N1y,    N1z = N1z                     , &
                                  N2x = N2x,    N2y = N2y,    N2z = N2z)
            END DO
         END DO
!$omp end parallel do

         CALL FASFMGCycle(this % Child,t,tol,lvl-1)
      END IF
!
!     ----------------------
!     Perform a V-Cycle here
!     ----------------------
!
      counter = 0
      IF (lvl > 1 ) THEN
         DO
            counter = counter + 1
            CALL FASVCycle(this,t,lvl,lvl)
            maxResidual = ComputeMaxResidual(this % p_sem)
            IF (maxval(maxResidual) <= tol) EXIT
         END DO
      ELSE
         DO
            counter = counter + 1
            dt = MaxTimeStep(this % p_sem, cfl )
            CALL SmoothIt   (this % p_sem, t, dt )
            maxResidual = ComputeMaxResidual(this % p_sem)
            
            IF (MOD(counter,100)==0) CALL PlotResiduals( lvl ,counter, this % p_sem)
            IF (maxval(maxResidual) <= tol) EXIT
         END DO
      END IF
      CALL PlotResiduals( lvl ,counter, this % p_sem ,.TRUE.)
!
!     --------------------------------------------------
!     If not on finest, Interpolate to next (finer) grid
!     --------------------------------------------------
! 
      IF (lvl < MGlevels) THEN
!$omp parallel do private(iEQ,N1x,N1y,N1z,N2x,N2y,N2z)
         DO iEl = 1, nelem
            DO iEQ = 1, N_EQN
               N1x = this % p_sem % Nx(iEl)
               N1y = this % p_sem % Ny(iEl)
               N1z = this % p_sem % Nz(iEl)
               N2x = this % Parent % p_sem % Nx(iEl)
               N2y = this % Parent % p_sem % Ny(iEl)
               N2z = this % Parent % p_sem % Nz(iEl)
               CALL Interpolate3D(Q1 = this          % p_sem % mesh % elements(iEl) % storage % Q(:,:,:,iEQ)      , &
                                  Q2 = this % Parent % p_sem % mesh % elements(iEl) % storage % Q(:,:,:,iEQ)      , &
                                  Interp = this % Prolongation(N2x,N2y,N2z) % Mat, &
                                  N1x = N1x,    N1y = N1y,    N1z = N1z             , &
                                  N2x = N2x,    N2y = N2y,    N2z = N2z)
            END DO
         END DO
!$omp end parallel do
      END IF
   END SUBROUTINE FASFMGCycle

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------
!  Subroutine that calls rhe smoother
!  ------------------------------------------
   SUBROUTINE MGRestrictToChild(this,lvl,t)
      IMPLICIT NONE
      !-------------------------------------------------------------
      CLASS(FASMultigrid_t), INTENT(INOUT) :: this     !<  Current level solver
      INTEGER              , INTENT(IN)    :: lvl
      REAL(KIND=RP)        , INTENT(IN)    :: t
      !-------------------------------------------------------------
      CLASS(FASMultigrid_t), POINTER       :: Child_p  ! The child
      INTEGER  :: iEl
      INTEGER  :: iEQ
      INTEGER  :: N1x,N1y,N1z
      INTEGER  :: N2x,N2y,N2z
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
               CALL Interpolate3D(Q1 = this    % p_sem % mesh % elements(iEl) % storage % Q(:,:,:,iEQ), &
                                  Q2 = Child_p % p_sem % mesh % elements(iEl) % storage % Q(:,:,:,iEQ), &
                                  Interp = this % Restriction(N1x,N1y,N1z) % Mat            , &
                                  N1x = N1x,    N1y = N1y,    N1z = N1z                     , &
                                  N2x = N2x,    N2y = N2y,    N2z = N2z)
            END DO
         END DO
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
               CALL Interpolate3D(Q1 = this    % p_sem % mesh % elements(iEl) % storage % Qdot(:,:,:,iEQ), &
                                  Q2 = Child_p % MGStorage(iEl) % S   (:,:,:,iEQ), &
                                  Interp = this % Restriction(N1x,N1y,N1z) % Mat    , &
                                  N1x = N1x,    N1y = N1y,    N1z = N1z             , &
                                  N2x = N2x,    N2y = N2y,    N2z = N2z)
            END DO
         END DO
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
         Child_p % tempsem % mesh % elements(iEl) % storage %S = 0._RP
      END DO
!$omp end parallel do
!
!     -------------------------------------------
!     If not on finest level, correct source term
!     -------------------------------------------
!      
      CALL ComputeTimeDerivative(Child_p % p_sem,t) 
      
!$omp parallel do
      DO iEl = 1, nelem
         Child_p % p_sem % mesh % elements(iEl) % storage % S = Child_p % MGStorage(iEl) % S - &
                                                                Child_p % p_sem % mesh % elements(iEl) % storage % Qdot
      END DO
!$omp end parallel do
      
      
   END SUBROUTINE MGRestrictToChild
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------
!  Subroutine that calls the smoother
!  ------------------------------------------
   SUBROUTINE SmoothIt(sem, t, dt)
      IMPLICIT NONE
      !---------------------------------------
      TYPE(DGSem)   :: sem
      REAL(KIND=RP) :: t, dt
      !---------------------------------------
      
      SELECT CASE (TRIM(SmootherName))
         CASE('RK3')
            CALL TakeRK3Step(sem, t, dt)
         CASE('SIRK')
!~            CALL TakeSIRK(sem, t, dt)
         CASE DEFAULT
            ERROR STOP 'Smoother not recognized'
      END SELECT
      
   END SUBROUTINE SmoothIt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------
!  Subroutine that destructs a FAS integrator
!  ------------------------------------------
   SUBROUTINE destruct(this)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(FASMultigrid_t), INTENT(INOUT) :: this
      !-----------------------------------------------------------
      
      CALL RecursiveDestructor(this,MGlevels)
      
   END SUBROUTINE destruct
   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   RECURSIVE SUBROUTINE RecursiveDestructor(Solver,lvl)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(FASMultigrid_t), INTENT(INOUT) :: Solver
      INTEGER                              :: lvl
      !-----------------------------------------------------------
      
      ! First go to finest level
      IF (lvl > 1) CALL RecursiveDestructor(Solver % Child,lvl-1)
      
      !Destruct Multigrid storage
      DEALLOCATE (Solver % MGStorage) ! allocatable components are automatically deallocated
      
      
      IF (lvl < MGlevels) THEN
         CALL Solver % tempsem % destruct()
         
         CALL Solver % p_sem % destruct()
         DEALLOCATE (Solver % p_sem)
         
         DEALLOCATE (Solver % Prolongation)
         
         NULLIFY    (Solver % Parent)
      ELSE
         NULLIFY    (Solver % p_sem)
      END IF
      
      IF (lvl > 1) THEN
         
         DEALLOCATE (Solver % Restriction)
         DEALLOCATE (Solver % Child)
      END IF
      
   END SUBROUTINE RecursiveDestructor
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
   SUBROUTINE CreateInterpolationOperators(Restriction,Prolongation,N1x,N1y,N1z,N2x,N2y,N2z,DeltaN,Lobatto)
      IMPLICIT NONE
      !-----------------------------------------------------
      TYPE(Interpolator_t), TARGET, intent(inout)  :: Restriction (0:,0:,0:)  !>  Restriction operator
      TYPE(Interpolator_t), TARGET, intent(inout)  :: Prolongation(0:,0:,0:)  !>  Prolongation operator
      INTEGER                     , intent(in)     :: N1x, N1y, N1z           !<  Fine grid order(anisotropic) of the element
      INTEGER                     , intent(out)    :: N2x, N2y, N2z           !>  Coarse grid order(anisotropic) of the element
      INTEGER                     , intent(in)     :: DeltaN                  !<  Interval of reduction of polynomial order for coarser level
      LOGICAL, OPTIONAL           , intent(in)     :: Lobatto                 !<  Is the quadrature a Legendre-Gauss-Lobatto representation?
      !-----------------------------------------------------
      TYPE(Interpolator_t), POINTER :: rest                    ! Pointer to constructed restriction interpolator
      TYPE(Interpolator_t), POINTER :: prol                    ! Pointer to constructed prolongation interpolator
      LOGICAL                       :: LGL = .FALSE.           ! Is the quadrature a Legendre-Gauss-Lobatto representation? (false is default)
      INTEGER                       :: i,j,k,l,m,n             ! Index counters
      INTEGER                       :: s,r                     ! Row/column counters for operators
      REAL(KIND=RP), ALLOCATABLE    :: x1 (:), y1 (:), z1 (:)  ! Position of quadrature points on mesh 1
      REAL(KIND=RP), ALLOCATABLE    :: w1x(:), w1y(:), w1z(:)  ! Weights for quadrature points on mesh 1
      REAL(KIND=RP), ALLOCATABLE    :: x2 (:), y2 (:), z2 (:)  ! Position of quadrature points on mesh 2
      REAL(KIND=RP), ALLOCATABLE    :: w2x(:), w2y(:), w2z(:)  ! Weights for quadrature points on mesh 2
      !-----------------------------------------------------
      
      IF (PRESENT(Lobatto) .AND. Lobatto) LGL = .TRUE.
!
!     --------------------------------------
!     Compute order of coarse representation
!     --------------------------------------
!
      N2x = N1x - DeltaN
      N2y = N1y - DeltaN
      N2z = N1z - DeltaN
      
      ! The order must be greater or equal to 0 (Legendre-Gauss quadrature) or 1 (Legendre-Gauss-Lobatto)
      IF (LGL) THEN
         IF (N2x < 1) N2x = 1
         IF (N2y < 1) N2y = 1
         IF (N2z < 1) N2z = 1
      ELSE
         IF (N2x < 1) N2x = 1       !! The threshold should actually be zero... Using 1 because max eigenvalue subroutine doesn't support 0
         IF (N2y < 1) N2y = 1
         IF (N2z < 1) N2z = 1
      END IF
      
      ! Return if the operators were already created
      IF (Restriction(N1x,N1y,N1z) % Created) RETURN
      
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
      IF (LGL) THEN
         CALL LegendreLobattoNodesAndWeights(N1x, x1, w1x)
         CALL LegendreLobattoNodesAndWeights(N1y, y1, w1y)
         CALL LegendreLobattoNodesAndWeights(N1z, z1, w1z)
         CALL LegendreLobattoNodesAndWeights(N2x, x2, w2x)
         CALL LegendreLobattoNodesAndWeights(N2y, y2, w2y)
         CALL LegendreLobattoNodesAndWeights(N2z, z2, w2z)
      ELSE
         CALL GaussLegendreNodesAndWeights(N1x, x1, w1x)
         CALL GaussLegendreNodesAndWeights(N1y, y1, w1y)
         CALL GaussLegendreNodesAndWeights(N1z, z1, w1z)
         CALL GaussLegendreNodesAndWeights(N2x, x2, w2x)
         CALL GaussLegendreNodesAndWeights(N2y, y2, w2y)
         CALL GaussLegendreNodesAndWeights(N2z, z2, w2z)
      END IF
!
!     -----------------------------
!     Fill the restriction operator
!     -----------------------------
!
      CALL Create3DRestrictionMatrix(rest % Mat,N1x,N1y,N1z,N2x,N2y,N2z,x1,y1,z1,x2,y2,z2,w1x,w1y,w1z,w2x,w2y,w2z)
!
!     ------------------------------
!     Fill the prolongation operator
!     ------------------------------
!
      CALL Create3DInterpolationMatrix(prol % Mat,N2x,N2y,N2z,N1x,N1y,N1z,x2,y2,z2,x1,y1,z1)
      
      ! All done
      rest % Created = .TRUE.
      prol % Created = .TRUE.
      
   END SUBROUTINE CreateInterpolationOperators
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------
!  Internal subroutine to print the residuals
!  ------------------------------------------
   subroutine PlotResiduals( iter, sweeps , sem, white )
      implicit none
      !--------------------------------------------------------
      integer    , intent(in)           :: iter
      TYPE(DGSem), intent(in)           :: sem
      integer    , intent(in)           :: sweeps
      LOGICAL    , intent(in), OPTIONAL :: white
      !--------------------------------------------------------
      real(kind=RP)             :: maxResiduals(N_EQN)
      CHARACTER(len=5)          :: color1
      CHARACTER(len=5)          :: color2
      !--------------------------------------------------------
      
      IF (PRESENT(white) .AND. white) THEN
         color1 = achar(27)//'[00m'
      else
         color1 = achar(27)//'[34m'
      END IF
      color2 = achar(27)//'[00m'
      
      IF( (MOD( ThisTimeStep+1, plotInterval) == 0) .or. (ThisTimeStep .eq. 0) ) THEN
         maxResiduals = ComputeMaxResidual(sem)
         write(STD_OUT , 110) color1,'FAS lvl/sweeps', iter ,"|",sweeps,"|", maxResiduals(IRHO) , "|" , maxResiduals(IRHOU) , &
                                 "|", maxResiduals(IRHOV) , "|" , maxResiduals(IRHOW) , "|" , maxResiduals(IRHOE),color2
      END IF
      
      110 format (A,A,I3,X,A,X,I10,X,A,X,ES10.3,X,A,X,ES10.3,X,A,X,ES10.3,X,A,X,ES10.3,X,A,X,ES10.3,A)
      
   end subroutine PlotResiduals
!
!/////////////////////////////////////////////////////////////////////////////////////////////////

END MODULE FASMultigridClass
