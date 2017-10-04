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
MODULE FASMultigridClass
   use FTValueDictionaryClass
   USE SMConstants
   USE ExplicitMethods
   USE DGSEMClass
   USE Physics
   
   USE PolynomialInterpAndDerivsModule
   USE GaussQuadrature
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
   
   ! Multigrid
   INTEGER        :: MGlevels       ! Total number of multigrid levels
   INTEGER        :: deltaN         ! 
   INTEGER        :: nelem          ! Number of elements (this is a p-multigrid implementation)
   INTEGER        :: ThisTimeStep   ! Current time step
   INTEGER        :: plotInterval   ! Read to display output
   LOGICAL        :: MGOutput       ! Display output?
   
   ! 
   LOGICAL        :: ManSol         ! Does this case have manufactured solutions?
   
   ! Variables related with snapshot files
   LOGICAL        :: WritePFiles = .FALSE.                   ! Write snapshot files?
   !! Manufactured solutions 11 elements
!~    INTEGER        :: PFileint = 20                          ! Iteration interval for writing snapshot files
!~    INTEGER        :: NPElems = 11                           ! Number of elements to track
!~    INTEGER        :: PElems(11) = &                         ! Elements to track
!~                      (/1,7,13,125,95,65,88,63,38,83,33/)
!~    INTEGER        :: PFileUnit(11)                     ! Units for writing files

!~    !! Manufactured solutions 3 elements
!~    INTEGER        :: PFileint = 60                          ! Iteration interval for writing snapshot files
!~    INTEGER        :: NPElems = 3                           ! Number of elements to track
!~    INTEGER        :: PElems(3) = &                         ! Elements to track
!~                      (/1,125,63/)
!~    INTEGER        :: PFileUnit(3)                     ! Units for writing files
   
   !! Flat plate
   INTEGER        :: PFileint = 50                        ! Iteration interval for writing snapshot files
   INTEGER        :: NPElems = 3                           ! Number of elements to track
   INTEGER        :: PElems(3) = &                         ! Elements to track
                     (/5,10,23/)
   INTEGER        :: PFileUnit(3) 
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
      
      IF (controlVariables % containsKey("cfl")) THEN
         cfl = controlVariables % doublePrecisionValueForKey("cfl")
      ELSE
         ERROR STOP '"cfl" keyword must be specified for the FAS integrator'
      END IF
      
!
!     -----------------------
!     Update module variables
!     -----------------------
!
      MGOutput       = controlVariables % logicalValueForKey("multigrid output")
      plotInterval   = controlVariables % integerValueForKey("output interval")
      ManSol         = sem % ManufacturedSol
      
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
      USE BoundaryConditionFunctions
      IMPLICIT NONE
      TYPE(FASMultigrid_t), TARGET  :: Solver
      INTEGER, DIMENSION(:)            :: N1x,N1y,N1z      !<  Order of approximation for every element in current solver
      INTEGER                          :: lvl              !<  Current multigrid level
      TYPE(FTValueDictionary)          :: controlVariables !< Control variables (for the construction of coarse sems
      !----------------------------------------------
      INTEGER, DIMENSION(nelem) :: N2x,N2y,N2z            !   Order of approximation for every element in child solver
      INTEGER                   :: N1xMAX,N1yMAX,N1zMAX   !   Maximum polynomial orders for current (fine) grid
      INTEGER                   :: N2xMAX,N2yMAX,N2zMAX   !   Maximum polynomial orders for child (coarse) grid
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
         Q1 = SIZE(Solver % p_sem % mesh % elements(k) % Q,1) - 1 
         Q2 = SIZE(Solver % p_sem % mesh % elements(k) % Q,2) - 1
         Q3 = SIZE(Solver % p_sem % mesh % elements(k) % Q,3) - 1
         Q4 = SIZE(Solver % p_sem % mesh % elements(k) % Q,4)
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
!~             IF(iEl==5) THEN
!~                print*, Solver % p_sem % Nx(iEl)
!~                print*, Solver % p_sem % Ny(iEl)
!~                print*, Solver % p_sem % Nz(iEl)
!~                print*,  Solver % MGStorage(iEl) % Scase
!~                read(*,*)
!~             ENDIF
         END DO
      END IF
      
!~       !/////////////////////////////////////////////////////////////////
!~       ! Temp block for plotting source term
!~       IF(lvl==2) THEN
!~       OPEN(newunit = fd, FILE = 'SourceInLowest.tec')
      
!~       WRITE(fd,'(A)') ' TITLE = "Free-stream preservation Euler equations" '
!~       WRITE(fd,'(A)') ' VARIABLES = "x","y","z","S_rho","S_rhou","S_rhov","S_rhow","S_rhoe"'
!~       DO iEl = 1, nelem
!~          Nxyz = Solver % p_sem % mesh % elements(iEl) % Nxyz
!~          WRITE(fd,*) "ZONE I=", Nxyz(1)+1, ",J=",Nxyz(2)+1, ",K=",Nxyz(3)+1,", F=POINT"
!~          DO k = 0, Nxyz(3)
!~             DO j = 0, Nxyz(2)
!~                DO i = 0, Nxyz(1)
                  
!~                   WRITE(fd,'(8E13.5)') Solver % p_sem % mesh % elements(iEl) % geom % x(1,i,j,k), &
!~                                        Solver % p_sem % mesh % elements(iEl) % geom % x(2,i,j,k), &
!~                                        Solver % p_sem % mesh % elements(iEl) % geom % x(3,i,j,k), &
!~                                        (Solver % MGStorage(iEl) % Scase (i,j,k,l), l = 1, 5)
!~                END DO
!~             END DO
!~          END DO 
!~       END DO
!~       STOP
!~       END IF
!~       !////////////////////////////////////////////////////////////////
      
      
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
         
         N2xMAX = N1xMAX - deltaN
         N2yMAX = N1yMAX - deltaN
         N2zMAX = N1zMAX - deltaN
         IF (N2xMAX < 1) N2xMAX = 1             ! TODO: Complete for Lobatto quadrature (max order = 1)
         IF (N2yMAX < 1) N2yMAX = 1             ! It is currently like it should be for LGL... Max eigenvalue routine!!!!
         IF (N2zMAX < 1) N2zMAX = 1
         
         ALLOCATE (Solver  % Restriction (0:N1xMAX,0:N1yMAX,0:N1zMAX))
         ALLOCATE (Child_p % Prolongation(0:N2xMAX,0:N2yMAX,0:N2zMAX))
         
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
         !Child_p % p_sem % ManufacturedSol = Solver % p_sem % ManufacturedSol
         
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
   SUBROUTINE solve(this,timestep,t)
      IMPLICIT NONE
      CLASS(FASMultigrid_t), INTENT(INOUT) :: this
      INTEGER                              :: timestep
      REAL(KIND=RP)                        :: t
      !-------------------------------------------------
      INTEGER                                 :: niter
      INTEGER                                 :: i
      CHARACTER(LEN=LINE_LENGTH)              :: FileName
      !-------------------------------------------------
      
      ThisTimeStep = timestep
!
!     --
!     Open pfiles for writing... And interpolate finest solution
!     --
!
      IF (WritePFiles) THEN
         IF (MOD(timestep+1,PFileint) == 0 .OR. timestep == 0) THEN
            
            DO i = 1, NPElems
               WRITE(FileName,'(A,I3.3,A,I4.4,A)') 'PFiles/Elem_', i ,'_' , timestep + 1,'.dat'
               OPEN( newunit = PFileUnit(i), FILE    = TRIM(FileName), ACTION  = 'WRITE')     
               
               WRITE(PFileUnit(i),*) timestep + 1                                                     ! Iteration
               WRITE(PFileUnit(i),*) this % p_sem % maxResidual                                       ! Gloval residual
               WRITE(PFileUnit(i),*) PElems(i)                                                        ! Element number
               WRITE(PFileUnit(i),*) MAXVAL(ABS(this % p_sem % mesh % elements(PElems(i)) % QDot))    ! Local residual
               
            END DO
            
            CALL InterpolateFinest(this,MGlevels)
         END IF
      END IF
!
!     -- 
!     Perform v-cycle
!     --
!
      
      CALL FASVCycle(this,t,MGlevels)
!
!     --
!     close opened files
!     --
!
      IF (WritePFiles) THEN
         IF (MOD(timestep+1,PFileint) == 0 .OR. timestep == 0) THEN
            DO i = 1, NPElems
               CLOSE( PFileUnit(i))   
            END DO
         END IF
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
               CALL Interpolate3D(Q1 = this    % tempsem % mesh % elements(iEl) % Q(:,:,:,iEQ), &
                                  Q2 = Child_p % tempsem % mesh % elements(iEl) % Q(:,:,:,iEQ), &
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
   RECURSIVE SUBROUTINE FASVCycle(this,t,lvl)
      IMPLICIT NONE
      !----------------------------------------------------------------------------
      CLASS(FASMultigrid_t), INTENT(INOUT) :: this    !<  Current level solver
      REAL(KIND=RP)                        :: t       !<  Simulation time
      INTEGER                              :: lvl     !<  Current multigrid level
      !----------------------------------------------------------------------------
      INTEGER                       :: iEl,iEQ        !Element/equation counter
      TYPE(FASMultigrid_t), POINTER :: Child_p        !Pointer to child
      INTEGER                       :: N1x, N1y, N1z  !Polynomial orders
      INTEGER                       :: N2x, N2y, N2z  !Polynomial orders
      REAL(KIND=RP)                 :: dt             !Time variables
      REAL(KIND=RP)                 :: maxResidual(N_EQN)   ! TODO: remove this from here... not needed
      !----------------------------------------------------------------------------
      ! Variables for writing PFiles
      INTEGER       :: ElemNmbr
      REAL(KIND=RP) :: TauNext, TauFinest
      !----------------------------------------------------------------------------
!
!     ------------------------------------
!     Copy fine grid solution to MGStorage
!        ... and clear source term
!     ------------------------------------
!
      IF (lvl < MGlevels) THEN
!$omp parallel do private(iEQ)
         DO iEl = 1, nelem
            this % MGStorage(iEl) % Q = this % p_sem % mesh % elements(iEl) % Q
            this % p_sem % mesh % elements(iEl) % S = 0._RP
            this % tempsem % mesh % elements(iEl) % S = 0._RP
         END DO
!$omp end parallel do
      END IF
!
!     -------------------------------------------
!     If not on finest level, correct source term
!     -------------------------------------------
!      
      IF (lvl < MGlevels) THEN
         CALL ComputeTimeDerivative(this % p_sem,t)
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Block that writes PFiles
         IF (WritePFiles) THEN
            IF (MOD(ThisTimeStep+1,PFileint) == 0 .OR. ThisTimeStep == 0) THEN
               
               CALL ComputeTimeDerivative(this % tempsem,t)
               
               DO iEl = 1, NPElems
                  ElemNmbr = PElems(iEl)
                  N1x = this % p_sem % Nx(ElemNmbr)
                  TauNext =   MAXVAL(ABS(this % p_sem   % mesh % elements(ElemNmbr) % Qdot + this % MGStorage(ElemNmbr) % Scase))
                  TauFinest = MAXVAL(ABS(this % tempsem % mesh % elements(ElemNmbr) % Qdot + this % MGStorage(ElemNmbr) % Scase))
                  WRITE(PFileUnit(iEl),*) lvl, N1x, TauFinest, TauNext
               END DO
            END IF
         END IF
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
!$omp parallel do
         DO iEl = 1, nelem
            this % p_sem % mesh % elements(iEl) % S = this % MGStorage(iEl) % S - &
                                                      this % p_sem % mesh % elements(iEl) % Qdot
         END DO
!$omp end parallel do
      END IF
!
!     -----------------------
!     Pre-smoothing procedure
!     -----------------------
!
      dt = MaxTimeStep(this % p_sem, cfl )
      CALL TakeRK3Step(this % p_sem, t, dt )
      IF (MGOutput) CALL PlotResiduals( lvl , this % p_sem )
      
      IF (lvl > 1) THEN
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
               CALL Interpolate3D(Q1 = this    % p_sem % mesh % elements(iEl) % Q(:,:,:,iEQ), &
                                  Q2 = Child_p % p_sem % mesh % elements(iEl) % Q(:,:,:,iEQ), &
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
               CALL Interpolate3D(Q1 = this    % p_sem % mesh % elements(iEl) % Qdot(:,:,:,iEQ), &
                                  Q2 = Child_p % MGStorage(iEl) % S   (:,:,:,iEQ), &
                                  Interp = this % Restriction(N1x,N1y,N1z) % Mat    , &
                                  N1x = N1x,    N1y = N1y,    N1z = N1z             , &
                                  N2x = N2x,    N2y = N2y,    N2z = N2z)
            END DO
         END DO
!$omp end parallel do
!
!        --------------------
!        Perform V-Cycle here
!        --------------------
!
         CALL FASVCycle(this % Child,t, lvl-1)
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
                                  Interp = Child_p % Prolongation(N1x,N1y,N1z) % Mat, &
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
            this % p_sem % mesh % elements(iEl) % Q = this % p_sem % mesh % elements(iEl) % Q + this % MGStorage(iEl) % E
         END DO
!$omp end parallel do
      
      END IF
!
!     ------------------------
!     Post-smoothing procedure
!     ------------------------
!
      dt = MaxTimeStep(this % p_sem, cfl )
      CALL TakeRK3Step(this % p_sem, t, dt)
      IF (MGOutput) CALL PlotResiduals( lvl , this % p_sem )
!
!     -------------------------
!     Compute coarse-grid error
!     -------------------------
!
      IF (lvl < MGlevels) THEN
!$omp parallel do private(iEQ)
         DO iEl = 1, nelem
            this % MGStorage(iEl) % E = this % p_sem % mesh % elements(iEl) % Q - this % MGStorage(iEl) % Q
         END DO
!$omp end parallel do
      END IF
      
   END SUBROUTINE FASVCycle
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --
!  TODO: finish this
!  --
   SUBROUTINE destruct(this)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(FASMultigrid_t), INTENT(INOUT) :: this
      !-----------------------------------------------------------
      
   END SUBROUTINE destruct
    
    
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
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
   SUBROUTINE CreateInterpolationOperators(Restriction,Prolongation,N1x,N1y,N1z,N2x,N2y,N2z,DeltaN,Lobatto)
      IMPLICIT NONE
!
!     ------------------------------------------------------------------------
!     Creates the restriction and prolongation operators for a certain element
!     for multigrid. Takes into account order anisotropy, but the coarse grid 
!     is constructed by reducing the polynomial order uniformly.
!     ------------------------------------------------------------------------
!
      !-----------------------------------------------------
      TYPE(Interpolator_t), TARGET  :: Restriction (0:,0:,0:)  !>  Restriction operator
      TYPE(Interpolator_t), TARGET  :: Prolongation(0:,0:,0:)  !>  Prolongation operator
      INTEGER                       :: N1x, N1y, N1z           !<  Fine grid order(anisotropic) of the element
      INTEGER                       :: N2x, N2y, N2z           !>  Coarse grid order(anisotropic) of the element
      INTEGER                       :: DeltaN                  !<  Interval of reduction of polynomial order for coarser level
      LOGICAL, OPTIONAL             :: Lobatto                 !<  Is the quadrature a Legendre-Gauss-Lobatto representation?
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
      prol => Prolongation(N2x,N2y,N2z)
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
!~       CALL Create3DInterpolationMatrix(rest % Mat,N1x,N1y,N1z,N2x,N2y,N2z,x1,y1,z1,x2,y2,z2)
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
!  ------
!  Internal subroutine to print the residuals
!  ----------------
   subroutine PlotResiduals( iter , sem )
      implicit none
      integer, intent(in)       :: iter
      TYPE(DGSem)               :: sem
      !--------------------------------------------------------
      real(kind=RP)             :: maxResiduals(N_EQN)
      !--------------------------------------------------------
      
      IF( (MOD( ThisTimeStep+1, plotInterval) == 0) .or. (ThisTimeStep .eq. 0) ) THEN
         maxResiduals = ComputeMaxResidual(sem)
         write(STD_OUT , 110) achar(27)//'[34mFAS level', iter ,"|            |", maxResiduals(IRHO) , "|" , maxResiduals(IRHOU) , &
                                 "|", maxResiduals(IRHOV) , "|" , maxResiduals(IRHOW) , "|" , maxResiduals(IRHOE),achar(27)//'[00m'
      END IF
      
      110 format (A,I8,X,A,X,ES10.3,X,A,X,ES10.3,X,A,X,ES10.3,X,A,X,ES10.3,X,A,X,ES10.3,A)
      
   end subroutine PlotResiduals
!
!/////////////////////////////////////////////////////////////////////////////////////////////////

END MODULE FASMultigridClass
