!! 
!////////////////////////////////////////////////////////////////////////
!
!      DGSEMClass.f95
!      Created: 2008-07-12 13:38:26 -0400 
!      By: David Kopriva
!
!      Basic class for the discontinuous Galerkin spectral element
!      solution of conservation laws.
!
!      Algorithms:
!         Algorithm 136: DGSEM Class
!         Algorithm 129: ConstructDGSem (Constructor)
!         Algorithm 138: ComputeTimeDerivative (TimeDerivative)
!         Algorithm 137: ComputeRiemannFluxes (EdgeFluxes)
!         Algorithm  35: InterpolateToFineMesh (2DCoarseToFineInterpolation)
!
!      Modified for 3D       6/11/15, 11:32 AM by DAK
!      Modified for mortars 25/04/17, 18:00    by arueda
!
!////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
Module DGSEMClass
   use SMConstants
   USE NodalStorageClass
   use MeshTypes
   use ElementClass
   USE HexMeshClass
   USE PhysicsStorage
#if defined(NAVIERSTOKES)
   USE ManufacturedSolutions
#endif
   use MonitorsClass
   use Physics
#ifdef _HAS_MPI_
   use mpi
#endif
   
   IMPLICIT NONE

   private
   public   ComputeQDot_FCN, DGSem, ConstructDGSem
   public   BCState_FCN, BCGradients_FCN

   public   DestructDGSEM, MaxTimeStep, ComputeMaxResiduals
   
   TYPE DGSem
      REAL(KIND=RP)                                           :: maxResidual
      integer                                                 :: nodes                 ! Either GAUSS or GAUSLOBATTO
      INTEGER                                                 :: numberOfTimeSteps
      INTEGER                                                 :: NDOF                         ! Number of degrees of freedom
      INTEGER           , ALLOCATABLE                         :: Nx(:), Ny(:), Nz(:)
      TYPE(HexMesh)                                           :: mesh
      PROCEDURE(BCState_FCN)    , NOPASS, POINTER :: externalState => NULL()
      PROCEDURE(BCGradients_FCN), NOPASS, POINTER :: externalGradients => NULL()
      LOGICAL                                                 :: ManufacturedSol = .FALSE.   ! Use manifactured solutions? default .FALSE.
      type(Monitor_t)                                         :: monitors
      contains
         procedure :: construct => ConstructDGSem
         procedure :: destruct  => DestructDGSem   
         procedure :: GetQ
         procedure :: SetQ
         procedure :: GetQdot
         procedure :: SaveSolutionForRestart
         procedure :: SetInitialCondition => DGSEM_SetInitialCondition
   END TYPE DGSem

   abstract interface
      SUBROUTINE BCState_FCN(x,t,nHat,Q,boundaryName)
         USE SMConstants
         REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
         REAL(KIND=RP)   , INTENT(INOUT) :: Q(:)
         CHARACTER(LEN=*), INTENT(IN)    :: boundaryName
      END SUBROUTINE BCState_FCN      

      SUBROUTINE BCGradients_FCN(x,t,nHat,gradU,boundaryName)
         USE SMConstants
         REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
         REAL(KIND=RP)   , INTENT(INOUT) :: gradU(:,:)
         CHARACTER(LEN=*), INTENT(IN)    :: boundaryName
      END SUBROUTINE BCGradients_FCN

      SUBROUTINE ComputeQDot_FCN( mesh, time, externalState, externalGradients )
         use SMConstants
         use HexMeshClass
         import BCState_FCN
         import BCGradients_FCN
         IMPLICIT NONE 
         type(HexMesh), target      :: mesh
         REAL(KIND=RP)              :: time
         procedure(BCState_FCN)     :: externalState
         procedure(BCGradients_FCN) :: externalGradients
      end subroutine ComputeQDot_FCN
   END INTERFACE

   CONTAINS 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructDGSem( self, meshFileName_, controlVariables, &
                                 externalState, externalGradients, polynomialOrder, Nx_, Ny_, Nz_, success, ChildSem )
      use ReadMeshFile
      use FTValueDictionaryClass
      use mainKeywordsModule
      use StopwatchClass
      use MPI_Process_Info
      use PartitionedMeshClass
      use MeshPartitioning
      IMPLICIT NONE
!
!     --------------------------
!     Constructor for the class.
!     --------------------------
!
      CLASS(DGSem)                       :: self                               !<> Class to be constructed
      character(len=*),         optional :: meshFileName_
      class(FTValueDictionary)           :: controlVariables                   !<  Name of mesh file
      EXTERNAL                           :: externalState, externalGradients   !<  External procedures that define the BCs
      INTEGER, OPTIONAL                  :: polynomialOrder(3)                 !<  Uniform polynomial order
      INTEGER, OPTIONAL, TARGET          :: Nx_(:), Ny_(:), Nz_(:)             !<  Non-uniform polynomial order
      LOGICAL, OPTIONAL                  :: success                            !>  Construction finalized correctly?
      logical, optional                  :: ChildSem                           !<  Is this a (multigrid) child sem?
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                     :: i,j,k,el                           ! Counters
      INTEGER, POINTER            :: Nx(:), Ny(:), Nz(:)                ! Orders of every element in mesh (used as pointer to use less space)
      integer                     :: nodes, NelL(2), NelR(2)
      INTEGER                     :: nTotalElem                              ! Number of elements in mesh
      INTEGER                     :: fUnit
      integer                     :: dir2D
      logical                     :: MeshInnerCurves                    ! The inner survaces of the mesh have curves?
      character(len=*), parameter :: TWOD_OFFSET_DIR_KEY = "2d mesh offset direction"
      INTERFACE
         SUBROUTINE externalState(x,t,nHat,Q,boundaryName)
            USE SMConstants
            REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
            REAL(KIND=RP)   , INTENT(INOUT) :: Q(:)
            CHARACTER(LEN=*), INTENT(IN)    :: boundaryName
         END SUBROUTINE externalState
         
         SUBROUTINE externalGradients(x,t,nHat,gradU,boundaryName)
            USE SMConstants
            REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
            REAL(KIND=RP)   , INTENT(INOUT) :: gradU(:,:)
            CHARACTER(LEN=*), INTENT(IN)    :: boundaryName
         END SUBROUTINE externalGradients
      END INTERFACE
      
      if ( present(ChildSem) .and. ChildSem ) self % mesh % child = .TRUE.
      
!
!     Measure preprocessing time
!     --------------------------      
      if (.not. self % mesh % child) then
         call Stopwatch % CreateNewEvent("Preprocessing")
         call Stopwatch % Start("Preprocessing")
      end if
      
      if ( present( meshFileName_ ) ) then
!
!        Mesh file set up by input argument
!        ----------------------------------
         self % mesh % meshFileName = trim(meshFileName_)

      else
!
!        Mesh file set up by controlVariables
!        ------------------------------------
         self % mesh % meshFileName = controlVariables % stringValueForKey(meshFileNameKey, requestedLength = LINE_LENGTH) 

      end if
!
!     ------------------------------
!     Discretization nodes selection
!     ------------------------------
!
      select case ( trim(controlVariables % stringValueForKey(trim(discretizationNodesKey), requestedLength = LINE_LENGTH)) )
      case("Gauss")
         nodes = GAUSS
      case("Gauss-Lobatto")
         nodes = GAUSSLOBATTO
      case default
         print*, "Unknown discretization nodes."
         print*, "Options available are:"
         print*, "   * Gauss"
         print*, "   * Gauss-Lobatto"
         errorMessage(STD_OUT)
         stop
      end select
      self % nodes = nodes
!
!     ---------------------------------------
!     Get polynomial orders for every element
!     ---------------------------------------
!
      IF (PRESENT(Nx_) .AND. PRESENT(Ny_) .AND. PRESENT(Nz_)) THEN
         Nx => Nx_
         Ny => Ny_
         Nz => Nz_
         nTotalElem = SIZE(Nx)
      ELSEIF (PRESENT(polynomialOrder)) THEN
         nTotalElem = NumOfElemsFromMeshFile( self % mesh % meshfileName )
         
         ALLOCATE (Nx(nTotalElem),Ny(nTotalElem),Nz(nTotalElem))
         Nx = polynomialOrder(1)
         Ny = polynomialOrder(2)
         Nz = polynomialOrder(3)
      ELSE
         ERROR STOP 'ConstructDGSEM: Polynomial order not specified'
      END IF
      
      ! Now store everything in sem
      IF (ALLOCATED(self % Nx)) DEALLOCATE (self % Nx)
      IF (ALLOCATED(self % Ny)) DEALLOCATE (self % Ny)
      IF (ALLOCATED(self % Nz)) DEALLOCATE (self % Nz)
      ALLOCATE (self % Nx(nTotalElem),self % Ny(nTotalElem),self % Nz(nTotalElem))
      self % Nx = Nx
      self % Ny = Ny
      self % Nz = Nz
      
      if ( max(maxval(Nx),maxval(Ny),maxval(Nz)) /= min(minval(Nx),minval(Ny),minval(Nz)) ) self % mesh % anisotropic = .TRUE.
      
!
!     -------------------------------------------------------------
!     Construct the polynomial storage for the elements in the mesh
!     -------------------------------------------------------------
!
      call NodalStorage(0) % Construct(nodes, 0)   ! Always construct orders 0 
      call NodalStorage(1) % Construct(nodes, 1)   ! and 1

      self % NDOF = 0
      DO k=1, nTotalElem
         self % NDOF = self % NDOF + N_EQN * (Nx(k) + 1) * (Ny(k) + 1) * (Nz(k) + 1)
         
         call NodalStorage(Nx(k)) % construct( nodes, Nx(k) )
         call NodalStorage(Ny(k)) % construct( nodes, Ny(k) )
         call NodalStorage(Nz(k)) % construct( nodes, Nz(k) )
      END DO
!
!     ------------------
!     Construct the mesh
!     ------------------
!
      if ( controlVariables % containsKey(TWOD_OFFSET_DIR_KEY) ) then
         select case ( controlVariables % stringValueForKey(TWOD_OFFSET_DIR_KEY,1))
         case("x")
            dir2D = 1
         case("y")
            dir2D = 2
         case("z")
            dir2D = 3
         case default
            print*, "Unrecognized 2D mesh offset direction"
            stop
            errorMessage(STD_OUT)
         end select

      else
         dir2D = 0

      end if
      
      if ( self % mesh % anisotropic .and. (min(minval(Nx),minval(Ny),minval(Nz)) <=1) .and. dir2D == 0 ) ERROR STOP ':: 3D anisotropic mesh must have N>=2'

      if (controlVariables % containsKey("mesh inner curves")) then
         MeshInnerCurves = controlVariables % logicalValueForKey("mesh inner curves")
      else
         MeshInnerCurves = .true.
      end if
!
!     **********************************************************
!     *                  MPI PREPROCESSING                     *
!     **********************************************************
!
!     Initialization
!     --------------
      call Initialize_MPI_Partitions
!
!     Prepare the processes to receive the partitions
!     -----------------------------------------------
      if ( MPI_Process % doMPIAction ) then
         call RecvPartitionMPI()
      end if
!
!     Read the mesh by the root rank to perform the partitioning
!     ----------------------------------------------------------
      if ( MPI_Process % doMPIRootAction ) then
!
!        Construct the full mesh
!        -----------------------
         call constructMeshFromFile( self % mesh, self % mesh % meshFileName, nodes, Nx, Ny, Nz, MeshInnerCurves , dir2D, success )
!
!        Perform the partitioning
!        ------------------------
         call PerformMeshPartitioning(self % mesh, MPI_Process % nProcs, mpi_allPartitions)
!
!        Send the partitions
!        -------------------
         call SendPartitionsMPI()
!
!        Destruct the full mesh
!        ----------------------
         call self % mesh % Destruct()

      end if
!
!     **********************************************************
!     *              MESH CONSTRUCTION                         *
!     **********************************************************
!
      CALL constructMeshFromFile( self % mesh, self % mesh % meshFileName, nodes, Nx, Ny, Nz, MeshInnerCurves , dir2D, success )
      
      IF(.NOT. success) RETURN
!
!     ------------------------
!     Allocate and zero memory
!     ------------------------
!
      DO k = 1, SIZE(self % mesh % elements) 
         CALL allocateElementStorage( self = self % mesh % elements(k), &
                                      nEqn = N_EQN, &
                                  nGradEqn = N_GRAD_EQN, &
                          computeGradients = computeGradients)
      END DO
!
!     ----------------------------------------------------
!     Get manufactured solution source term (if requested)
!     ----------------------------------------------------
!
#if defined(NAVIERSTOKES)
      IF (self % ManufacturedSol) THEN
         DO el = 1, SIZE(self % mesh % elements) 
            DO k=0, Nz(el)
               DO j=0, Ny(el)
                  DO i=0, Nx(el)
                     IF (flowIsNavierStokes) THEN
                        CALL ManufacturedSolutionSourceNS(self % mesh % elements(el) % geom % x(:,i,j,k), &
                                                          0._RP, &
                                                          self % mesh % elements(el) % storage % S (:,i,j,k)  )
                     ELSE
                        CALL ManufacturedSolutionSourceEuler(self % mesh % elements(el) % geom % x(:,i,j,k), &
                                                             0._RP, &
                                                             self % mesh % elements(el) % storage % S (:,i,j,k)  )
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END IF
#endif
!
!     -----------------------
!     Set boundary conditions
!     -----------------------
!
      self % externalState     => externalState
      self % externalGradients => externalGradients
      
      call assignBoundaryConditions(self)
!
!     ------------------
!     Build the monitors
!     ------------------
!
      self % monitors = ConstructMonitors(self % mesh, controlVariables)
      
      NULLIFY(Nx,Ny,Nz)
!
!     Stop measuring preprocessing time
!     ----------------------------------
      if (.not. self % mesh % child) call Stopwatch % Pause("Preprocessing")

      END SUBROUTINE ConstructDGSem

!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructDGSem( self )
      IMPLICIT NONE 
      CLASS(DGSem) :: self
      INTEGER      :: k      !Counter
      
      CALL self % mesh % destruct
      self % externalState     => NULL()
      self % externalGradients => NULL()
      END SUBROUTINE DestructDGSem
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SaveSolutionForRestart( self, fUnit ) 
         IMPLICIT NONE
         CLASS(DGSem)     :: self
         INTEGER          :: fUnit
         INTEGER          :: k

         DO k = 1, SIZE(self % mesh % elements) 
            WRITE(fUnit) self % mesh % elements(k) % storage % Q
         END DO

      END SUBROUTINE SaveSolutionForRestart

      subroutine DGSEM_SetInitialCondition( self, controlVariables, initial_iteration, initial_time ) 
         use FTValueDictionaryClass
         USE mainKeywordsModule
         implicit none
         class(DGSEM)   :: self
         class(FTValueDictionary), intent(in)   :: controlVariables
         integer                                :: restartUnit
         integer,       intent(out)             :: initial_iteration
         real(kind=RP), intent(out)             :: initial_time 
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=LINE_LENGTH)             :: fileName, solutionName
         logical                                :: saveGradients
         interface
            SUBROUTINE UserDefinedInitialCondition(mesh, thermodynamics_, &
                                                           dimensionless_,&
                                                           refValues_)
               USE SMConstants
               use PhysicsStorage
               use HexMeshClass
               implicit none
               class(HexMesh)                  :: mesh
               type(Thermodynamics_t), intent(in)  :: thermodynamics_
               type(Dimensionless_t),  intent(in)  :: dimensionless_
               type(RefValues_t),      intent(in)  :: refValues_
            END SUBROUTINE UserDefinedInitialCondition
            character(len=LINE_LENGTH) function getFileName( inputLine )
               use SMConstants
               implicit none
               character(len=*)     :: inputLine
            end function getFileName
         end interface

         IF ( controlVariables % logicalValueForKey(restartKey) )     THEN
            fileName = controlVariables % stringValueForKey(restartFileNameKey,requestedLength = LINE_LENGTH)
            CALL self % mesh % LoadSolution(fileName, initial_iteration, initial_time)
         ELSE
   
            call UserDefinedInitialCondition(self % mesh, thermodynamics, &
                                                    dimensionless, &
                                                        refValues )
            initial_time = 0.0_RP
            initial_iteration = 0
!
!           Save the initial condition (for initialized cases)
!           --------------------------
            saveGradients = controlVariables % logicalValueForKey(saveGradientsToSolutionKey)
            solutionName = controlVariables % stringValueForKey(solutionFileNameKey, requestedLength = LINE_LENGTH)
            solutionName = trim(getFileName(solutionName))
            write(solutionName,'(A,A,I10.10,A)') trim(solutionName), "_", initial_iteration, ".hsol"
            call self % mesh % SaveSolution(initial_iteration, initial_time, solutionName, saveGradients)

         END IF
   
      end subroutine DGSEM_SetInitialCondition
!
!//////////////////////////////////////////////////////////////////////// 
!
!  Routine to set the solution in each element with a global solution vector
!
   SUBROUTINE SetQ(self,Q)
      IMPLICIT NONE
      CLASS(DGSem)   ,     INTENT(INOUT)           :: self 
      REAL(KIND = RP),     INTENT(IN)              :: Q(:)   
      
      INTEGER                                      :: Nx, Ny, Nz, l, i, j, k, counter, elm
      
      IF (SIZE(Q) /= self % NDOF) ERROR STOP 'Size mismatch in DGSEM:SetQ'
      
      counter = 1
      DO elm = 1, size(self%mesh%elements)
         Nx = self%mesh%elements(elm)%Nxyz(1)
         Ny = self%mesh%elements(elm)%Nxyz(2)
         Nz = self%mesh%elements(elm)%Nxyz(3)
         DO k = 0, Nz
            DO j = 0, Ny
               DO i = 0, Nx
                  DO l = 1,N_EQN
                     self%mesh%elements(elm)%storage%Q(l,i,j,k) = Q(counter) ! This creates a temporary array: storage must be modified to avoid that
                     counter =  counter + 1
                  END DO
               END DO
            END DO
         END DO
      END DO 
         
   END SUBROUTINE SetQ
!
!////////////////////////////////////////////////////////////////////////////////////////       
!
!  Routine to get the solution in each element as a global solution vector
!
   SUBROUTINE GetQ(self,Q)
      IMPLICIT NONE
      CLASS(DGSem),        INTENT(INOUT)            :: self
      REAL(KIND = RP),     INTENT(OUT)              :: Q(:)
      
      INTEGER                                       :: Nx, Ny, Nz, l, i, j, k, counter, elm
      
      IF (SIZE(Q) /= self % NDOF) ERROR STOP 'Size mismatch in DGSEM:GetQ'
      counter = 1
      DO elm = 1, size(self%mesh%elements)
         Nx = self%mesh%elements(elm)%Nxyz(1)
         Ny = self%mesh%elements(elm)%Nxyz(2)
         Nz = self%mesh%elements(elm)%Nxyz(3)
         DO k = 0, Nz
            DO j = 0, Ny
                DO i = 0, Nx
                  DO l = 1,N_EQN
                     Q(counter)  = self%mesh%elements(elm)%storage%Q(l,i, j, k) ! This creates a temporary array: storage must be modified to avoid that
                     counter =  counter + 1
                  END DO
                END DO
            END DO
         END DO
      END DO
      
   END SUBROUTINE GetQ
!
!////////////////////////////////////////////////////////////////////////////////////////      
!
!  Routine to get the solution's time derivative in each element as a global solution vector
!
   SUBROUTINE GetQdot(self,Qdot)
      IMPLICIT NONE
      CLASS(DGSem),        INTENT(INOUT)            :: self
      REAL(KIND = RP),     INTENT(OUT)              :: Qdot(:)
      
      INTEGER                                       :: Nx, Ny, Nz, l, i, j, k, counter, elm
      
      IF (SIZE(Qdot) /= self % NDOF) ERROR STOP 'Size mismatch in DGSEM:GetQdot'
      counter = 1
      DO elm = 1, size(self%mesh%elements)
         Nx = self%mesh%elements(elm)%Nxyz(1)
         Ny = self%mesh%elements(elm)%Nxyz(2)
         Nz = self%mesh%elements(elm)%Nxyz(3)
         DO k = 0, Nz
            DO j = 0, Ny
               DO i = 0, Nx
                  DO l = 1,N_EQN
                     Qdot(counter)  = self%mesh%elements(elm)%storage%Qdot(l,i, j, k) ! This creates a temporary array: storage must be modified to avoid that
                     counter =  counter + 1
                  END DO
               END DO
            END DO
         END DO
      END DO
      
   END SUBROUTINE GetQdot
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE assignBoundaryConditions(self)
!
!        ------------------------------------------------------------
!        Assign the boundary condition type to the boundaries through
!        their boundary names
!        ------------------------------------------------------------
!
         USE SharedBCModule
         IMPLICIT NONE  
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(DGSem)     :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                         :: eID, k
         CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryType, boundaryName
         
         DO eID = 1, SIZE( self % mesh % elements)
            DO k = 1,6
               boundaryName = self % mesh % elements(eID) % boundaryName(k)
               IF ( boundaryName /= emptyBCName )     THEN
                  boundaryType = bcTypeDictionary % stringValueForKey(key             = boundaryName, &
                                                                      requestedLength = BC_STRING_LENGTH)
                  IF( LEN_TRIM(boundaryType) > 0) then
                     self % mesh % elements(eID) % boundaryType(k) = boundaryType ! TODO bType and bName in elements are needed?
                     self % mesh % faces(self % mesh % elements(eID) % faceIDs(k)) % boundaryType = boundaryType
                  end if
               END IF 
                                                                   
            END DO  
         END DO

      END SUBROUTINE assignBoundaryConditions
!
!////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------
!  Compute maximum residual L_inf norm
!  -----------------------------------
!
#if defined(NAVIERSTOKES)
   FUNCTION ComputeMaxResiduals(mesh) RESULT(maxResidual)
      use MPI_Process_Info
      IMPLICIT NONE
      !----------------------------------------------
      CLASS(HexMesh), intent(in)  :: mesh
      REAL(KIND=RP) :: maxResidual(N_EQN)
      !----------------------------------------------
      INTEGER       :: id , eq, ierr
      REAL(KIND=RP) :: localMaxResidual(N_EQN)
      real(kind=RP) :: localRho, localRhou, localRhov, localRhow, localRhoe
      real(kind=RP) :: Rho, Rhou, Rhov, Rhow, Rhoe
      !----------------------------------------------
      
      maxResidual = 0.0_RP
      Rho = 0.0_RP
      Rhou = 0.0_RP
      Rhov = 0.0_RP
      Rhow = 0.0_RP
      Rhoe = 0.0_RP

!$omp parallel shared(maxResidual, Rho, Rhou, Rhov, Rhow, Rhoe, mesh) default(private)
!$omp do reduction(max:Rho,Rhou,Rhov,Rhow,Rhoe)
      DO id = 1, SIZE( mesh % elements )
         localRho = maxval(abs(mesh % elements(id) % storage % QDot(IRHO,:,:,:)))
         localRhou = maxval(abs(mesh % elements(id) % storage % QDot(IRHOU,:,:,:)))
         localRhov = maxval(abs(mesh % elements(id) % storage % QDot(IRHOV,:,:,:)))
         localRhow = maxval(abs(mesh % elements(id) % storage % QDot(IRHOW,:,:,:)))
         localRhoe = maxval(abs(mesh % elements(id) % storage % QDot(IRHOE,:,:,:)))
      
         Rho = max(Rho,localRho)
         Rhou = max(Rhou,localRhou)
         Rhov = max(Rhov,localRhov)
         Rhow = max(Rhow,localRhow)
         Rhoe = max(Rhoe,localRhoe)
      END DO
!$omp end do
!$omp end parallel

      maxResidual = (/Rho, Rhou, Rhov, Rhow, Rhoe/)

#ifdef _HAS_MPI_
      if ( MPI_Process % doMPIAction ) then
         localMaxResidual = maxResidual
         call mpi_allreduce(localMaxResidual, maxResidual, N_EQN, MPI_DOUBLE, MPI_MAX, &
                            MPI_COMM_WORLD, ierr)
      end if
#endif

   END FUNCTION ComputeMaxResiduals
#elif defined(CAHNHILLIARD)
   FUNCTION ComputeMaxResiduals(mesh) RESULT(maxResidual)
      use MPI_Process_Info
      IMPLICIT NONE
      !----------------------------------------------
      CLASS(HexMesh)  :: mesh
      REAL(KIND=RP) :: maxResidual(N_EQN)
      !----------------------------------------------
      INTEGER       :: id , eq, ierr
      real(kind=RP) :: cMax
      REAL(KIND=RP) :: localCMax(NCONS)
      !----------------------------------------------
      
      maxResidual = 0.0_RP

      cMax = 0.0_RP
      localCMax = 0.0_RP
!$omp parallel shared(maxResidual,cMax,mesh) default(private)
!$omp do reduction(max:cMax)
      DO id = 1, SIZE( mesh % elements )
         localCMax(1) = maxval(abs(mesh % elements(id) % storage % QDot(1,:,:,:)))
         
         cMax = max(localCMax(1), cMax)
      END DO
!$omp end do
!$omp end parallel

      maxResidual = cMax

#ifdef _HAS_MPI_
      if ( MPI_Process % doMPIAction ) then
         localCMax = maxResidual
         call mpi_allreduce(localCMax, maxResidual, N_EQN, MPI_DOUBLE, MPI_MAX, &
                            MPI_COMM_WORLD, ierr)
      end if
#endif
   END FUNCTION ComputeMaxResiduals
#endif
!
!//////////////////////////////////////////////////////////////////////// 
!
!  -------------------------------------------------------------------
!  Estimate the maximum time-step of the system. This
!  routine computes a heuristic based on the smallest mesh
!  spacing (which goes as 1/N^2) AND the eigenvalues of the
!  particular hyperbolic system being solved (Euler or Navier Stokes). This is not
!  the only way to estimate the time-step, but it works in practice.
!  Other people use variations on this and we make no assertions that
!  it is the only or best way. Other variations look only at the smallest
!  mesh values, other account for differences across the element.
!     -------------------------------------------------------------------
   real(kind=RP) function MaxTimeStep( self, cfl, dcfl )
      use VariableConversion
      use MPI_Process_Info
      implicit none
      !------------------------------------------------
      type(DGSem)                :: self
      real(kind=RP), intent(in)  :: cfl      !<  Advective cfl number
      real(kind=RP), optional, intent(in)  :: dcfl     !<  Diffusive cfl number
      !------------------------------------------------
      integer                       :: i, j, k, eID                     ! Coordinate and element counters
      integer                       :: N(3)                             ! Polynomial order in the three reference directions
      integer                       :: ierr                             ! Error for MPI calls
      real(kind=RP)                 :: eValues(3)                       ! Advective eigenvalues
      real(kind=RP)                 :: dcsi, deta, dzet                 ! Smallest mesh spacing
      real(kind=RP)                 :: dcsi2, deta2, dzet2              ! Smallest mesh spacing squared
      real(kind=RP)                 :: lamcsi_a, lamzet_a, lameta_a     ! Advective eigenvalues in the three reference directions
      real(kind=RP)                 :: lamcsi_v, lamzet_v, lameta_v     ! Diffusive eigenvalues in the three reference directions
      real(kind=RP)                 :: jac, mu, T                       ! Mapping Jacobian, viscosity and temperature
      real(kind=RP)                 :: Q(N_EQN)                         ! The solution in a node
      real(kind=RP)                 :: TimeStep_Conv, TimeStep_Visc     ! Time-step for convective and diffusive terms
      real(kind=RP)                 :: localMax_dt_v, localMax_dt_a     ! Time step to perform MPI reduction
      type(NodalStorage_t), pointer :: spAxi_p, spAeta_p, spAzeta_p     ! Pointers to the nodal storage in every direction
      external                      :: ComputeEigenvaluesForState       ! Advective eigenvalues
      !--------------------------------------------------------
#if defined(NAVIERSTOKES)      
!     Initializations
!     ---------------
      
      TimeStep_Conv = huge(1._RP)
      TimeStep_Visc = huge(1._RP)
      
!$omp parallel shared(self,TimeStep_Conv,TimeStep_Visc,NodalStorage,cfl,dcfl,flowIsNavierStokes) default(private) 
!$omp do reduction(min:TimeStep_Conv,TimeStep_Visc)
      do eID = 1, SIZE(self % mesh % elements) 
         N = self % mesh % elements(eID) % Nxyz
         spAxi_p => NodalStorage(N(1))
         spAeta_p => NodalStorage(N(2))
         spAzeta_p => NodalStorage(N(3))
         
         if ( N(1) .ne. 0 ) then
            dcsi = 1.0_RP / abs( spAxi_p   % x(1) - spAxi_p   % x (0) )   

         else
            dcsi = 0.0_RP

         end if

         if ( N(2) .ne. 0 ) then
            deta = 1.0_RP / abs( spAeta_p  % x(1) - spAeta_p  % x (0) )
         
         else
            deta = 0.0_RP

         end if

         if ( N(3) .ne. 0 ) then
            dzet = 1.0_RP / abs( spAzeta_p % x(1) - spAzeta_p % x (0) )

         else
            dzet = 0.0_RP

         end if
         
         if (flowIsNavierStokes) then
            dcsi2 = dcsi*dcsi
            deta2 = deta*deta
            dzet2 = dzet*dzet
         end if
         
         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
!
!           ------------------------------------------------------------
!           The maximum eigenvalues for a particular state is determined
!           by the physics.
!           ------------------------------------------------------------
!
            Q(1:N_EQN) = self % mesh % elements(eID) % storage % Q(1:N_EQN,i,j,k)
            CALL ComputeEigenvaluesForState( Q , eValues )
            
            jac      = self % mesh % elements(eID) % geom % jacobian(i,j,k)
!
!           ----------------------------
!           Compute contravariant values
!           ----------------------------
!              
            lamcsi_a =abs( self % mesh % elements(eID) % geom % jGradXi(IX,i,j,k)   * eValues(IX) + &
                           self % mesh % elements(eID) % geom % jGradXi(IY,i,j,k)   * eValues(IY) + &
                           self % mesh % elements(eID) % geom % jGradXi(IZ,i,j,k)   * eValues(IZ) ) * dcsi
  
            lameta_a =abs( self % mesh % elements(eID) % geom % jGradEta(IX,i,j,k)  * eValues(IX) + &
                           self % mesh % elements(eID) % geom % jGradEta(IY,i,j,k)  * eValues(IY) + &
                           self % mesh % elements(eID) % geom % jGradEta(IZ,i,j,k)  * eValues(IZ) ) * deta
  
            lamzet_a =abs( self % mesh % elements(eID) % geom % jGradZeta(IX,i,j,k) * eValues(IX) + &
                           self % mesh % elements(eID) % geom % jGradZeta(IY,i,j,k) * eValues(IY) + &
                           self % mesh % elements(eID) % geom % jGradZeta(IZ,i,j,k) * eValues(IZ) ) * dzet
            
            TimeStep_Conv = min( TimeStep_Conv, cfl*abs(jac)/(lamcsi_a+lameta_a+lamzet_a) )
            
            if (flowIsNavierStokes) then
               T        = Temperature(Q)
               mu       = SutherlandsLaw(T)
               lamcsi_v = mu * dcsi2 * abs(sum(self % mesh % elements(eID) % geom % jGradXi  (:,i,j,k)))
               lameta_v = mu * deta2 * abs(sum(self % mesh % elements(eID) % geom % jGradEta (:,i,j,k)))
               lamzet_v = mu * dzet2 * abs(sum(self % mesh % elements(eID) % geom % jGradZeta(:,i,j,k)))
               
               TimeStep_Visc = min( TimeStep_Visc, dcfl*abs(jac)/(lamcsi_v+lameta_v+lamzet_v) )
            end if
                  
         end do ; end do ; end do
      end do 
!$omp end do
!$omp end parallel
         
#ifdef _HAS_MPI_
      if ( MPI_Process % doMPIAction ) then
         localMax_dt_v = TimeStep_Visc
         localMax_dt_a = TimeStep_Conv
         call mpi_allreduce(localMax_dt_a, TimeStep_Conv, 1, MPI_DOUBLE, MPI_MIN, &
                            MPI_COMM_WORLD, ierr)
         call mpi_allreduce(localMax_dt_v, TimeStep_Visc, 1, MPI_DOUBLE, MPI_MIN, &
                            MPI_COMM_WORLD, ierr)
      end if
#endif
      
      if (TimeStep_Conv  < TimeStep_Visc) then
         self % mesh % dt_restriction = DT_CONV
         MaxTimeStep  = TimeStep_Conv
      else
         self % mesh % dt_restriction = DT_DIFF
         MaxTimeStep  = TimeStep_Visc
      end if
#endif
   end function MaxTimeStep
!
!////////////////////////////////////////////////////////////////////////
!
end module DGSEMClass
