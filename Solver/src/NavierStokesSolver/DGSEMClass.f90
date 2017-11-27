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
   USE NodalStorageClass
   USE HexMeshClass
   USE PhysicsStorage
   USE SpatialDiscretization
   USE ManufacturedSolutions
   use MonitorsClass
   
   IMPLICIT NONE
   
   ABSTRACT INTERFACE
      SUBROUTINE externalStateSubroutine(x,t,nHat,Q,boundaryName)
         USE SMConstants
         REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
         REAL(KIND=RP)   , INTENT(INOUT) :: Q(:)
         CHARACTER(LEN=*), INTENT(IN)    :: boundaryName
      END SUBROUTINE externalStateSubroutine
      
      SUBROUTINE externalGradientsSubroutine(x,t,nHat,gradU,boundaryName)
         USE SMConstants
         REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
         REAL(KIND=RP)   , INTENT(INOUT) :: gradU(:,:)
         CHARACTER(LEN=*), INTENT(IN)    :: boundaryName
      END SUBROUTINE externalGradientsSubroutine
   END INTERFACE
   
   TYPE DGSem
      REAL(KIND=RP)                                           :: maxResidual
      INTEGER                                                 :: numberOfTimeSteps
      INTEGER                                                 :: NDOF                         ! Number of degrees of freedom
      TYPE(NodalStorage)                        , allocatable :: spA(:)
      INTEGER           , ALLOCATABLE                         :: Nx(:), Ny(:), Nz(:)
      TYPE(HexMesh)                                           :: mesh
      PROCEDURE(externalStateSubroutine)    , NOPASS, POINTER :: externalState => NULL()
      PROCEDURE(externalGradientsSubroutine), NOPASS, POINTER :: externalGradients => NULL()
      LOGICAL                                                 :: ManufacturedSol = .FALSE.   ! Use manifactured solutions? default .FALSE.
      type(Monitor_t)                                        :: monitors
      contains
         procedure :: construct => ConstructDGSem
         procedure :: destruct  => DestructDGSem   
         procedure :: GetQ
         procedure :: SetQ
         procedure :: GetQdot
         procedure :: SaveSolutionForRestart
         procedure :: SetInitialCondition => DGSEM_SetInitialCondition
   END TYPE DGSem
   
      
   CONTAINS 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructDGSem( self, meshFileName_, controlVariables, &
                                 externalState, externalGradients, polynomialOrder, Nx_, Ny_, Nz_, success )
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
      character(len=LINE_LENGTH)  :: meshFileName
      logical                     :: MeshInnerCurves                    ! The inner survaces of the mesh have curves?
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
!
!     Measure preprocessing time
!     --------------------------      
      call Stopwatch % CreateNewEvent("Preprocessing")
      call Stopwatch % Start("Preprocessing")

      if ( present( meshFileName_ ) ) then
!
!        Mesh file set up by input argument
!        ----------------------------------
         meshFileName = trim(meshFileName_)

      else
!
!        Mesh file set up by controlVariables
!        ------------------------------------
         meshFileName = controlVariables % stringValueForKey(meshFileNameKey, requestedLength = LINE_LENGTH) 

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
         nTotalElem = NumOfElemsFromMeshFile( meshfileName )
         
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
!
!     -------------------------------------------------------------
!     Construct the polynomial storage for the elements in the mesh
!     -------------------------------------------------------------
!
      IF (ALLOCATED(self % spA)) DEALLOCATE(self % spa)
      k = max( MAXVAL(Nx), MAXVAL(Ny), MAXVAL(Nz) )
      ALLOCATE(self % spA(0:k) )
      
      self % NDOF = 0
      DO k=1, nTotalElem
         self % NDOF = self % NDOF + N_EQN * (Nx(k) + 1) * (Ny(k) + 1) * (Nz(k) + 1)
         
         call self % spA(Nx(k)) % construct( nodes, Nx(k) )
         call self % spA(Ny(k)) % construct( nodes, Ny(k) )
         call self % spA(Nz(k)) % construct( nodes, Nz(k) )
      END DO
!
!     ------------------
!     Construct the mesh
!     ------------------
!
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
         call constructMeshFromFile( self % mesh, meshfileName, nodes, self % spA, Nx, Ny, Nz, MeshInnerCurves , success )
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
      CALL constructMeshFromFile( self % mesh, meshfileName, nodes, self % spA, Nx, Ny, Nz, MeshInnerCurves , success )
      
      IF(.NOT. success) RETURN
!
!     ------------------------
!     Allocate and zero memory
!     ------------------------
!
      DO k = 1, SIZE(self % mesh % elements) 
         CALL allocateElementStorage( self % mesh % elements(k), Nx(k),Ny(k),Nz(k), &
                                      N_EQN, N_GRAD_EQN, flowIsNavierStokes )
      END DO
!
!     ----------------------------------------------------
!     Get manufactured solution source term (if requested)
!     ----------------------------------------------------
!
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
!
!     -----------------------------------------
!     Initialize Spatial discretization methods
!     -----------------------------------------
!
      call Initialize_SpaceAndTimeMethods(controlVariables)
      
      NULLIFY(Nx,Ny,Nz)
!
!     Stop measuring preprocessing time
!     ----------------------------------
      call Stopwatch % Pause("Preprocessing")

      END SUBROUTINE ConstructDGSem

!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructDGSem( self )
      IMPLICIT NONE 
      CLASS(DGSem) :: self
      INTEGER      :: k      !Counter
      
      DO k=0, UBOUND(self % spA,1)
         IF (.NOT. self % spA(k) % Constructed) CYCLE
         CALL self % spA(k) % destruct()
      END DO
      
      CALL self % mesh % destruct
      self % externalState     => NULL()
      self % externalGradients => NULL()
      IF ( ALLOCATED(InviscidMethod) ) DEALLOCATE( InviscidMethod )
      IF ( ALLOCATED(ViscousMethod ) ) DEALLOCATE( ViscousMethod ) 
      
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
!////////////////////////////////////////////////////////////////////////////////////////      
!
!  -----------------------------------
!  Compute maximum residual L_inf norm
!  -----------------------------------
   FUNCTION ComputeMaxResidual(self) RESULT(maxResidual)
      use MPI_Process_Info
#ifdef _HAS_MPI_
      use mpi
#endif
      IMPLICIT NONE
      !----------------------------------------------
      CLASS(DGSem)  :: self
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

!$omp parallel shared(maxResidual, Rho, Rhou, Rhov, Rhow, Rhoe, self) default(private)
!$omp do reduction(max:Rho,Rhou,Rhov,Rhow,Rhoe)
      DO id = 1, SIZE( self % mesh % elements )
         localRho = maxval(abs(self % mesh % elements(id) % storage % QDot(IRHO,:,:,:)))
         localRhou = maxval(abs(self % mesh % elements(id) % storage % QDot(IRHOU,:,:,:)))
         localRhov = maxval(abs(self % mesh % elements(id) % storage % QDot(IRHOV,:,:,:)))
         localRhow = maxval(abs(self % mesh % elements(id) % storage % QDot(IRHOW,:,:,:)))
         localRhoe = maxval(abs(self % mesh % elements(id) % storage % QDot(IRHOE,:,:,:)))
      
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
#endif _HAS_MPI_

   END FUNCTION ComputeMaxResidual
!
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
      SUBROUTINE ComputeTimeDerivative( self, time )
         USE SpatialDiscretization
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(DGSem)   :: self
         REAL(KIND=RP) :: time
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: k
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel shared(self, time)
         call self % mesh % ProlongSolutionToFaces(self % spA)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
         call self % mesh % UpdateMPIFacesSolution
!
!        -----------------
!        Compute gradients
!        -----------------
!
         if ( flowIsNavierStokes ) then
            CALL DGSpatial_ComputeGradient( self % mesh , self % spA, time , self % externalState , self % externalGradients )
         end if
!
!        -------------------------------------------------------
!        Inviscid Riemann fluxes from the solutions on the faces
!        -------------------------------------------------------
!
         CALL ComputeRiemannFluxes( self, time )
!
!        -----------------------
!        Compute time derivative
!        -----------------------
!
         call TimeDerivative_ComputeQDot( self % mesh , time )
!$omp end parallel
!
      END SUBROUTINE ComputeTimeDerivative
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE computeRiemannFluxes( self, time )
         USE Physics
         USE BoundaryConditionFunctions
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(DGSem)   :: self
         REAL(KIND=RP) :: time
!
!        ---------------
!        Local Variables
!        ---------------
!
         INTEGER       :: faceID
         INTEGER       :: eIDLeft, eIDRight
         INTEGER       :: fIDLeft

!$omp do schedule(runtime)
         DO faceID = 1, SIZE( self % mesh % faces)
            associate( f => self % mesh % faces(faceID))
            select case (f % faceType)
            case (HMESH_INTERIOR)
               CALL computeElementInterfaceFlux ( f )

            case (HMESH_BOUNDARY)
               CALL computeBoundaryFlux(f, time, &
                                        self % externalState , self % externalGradients)

            case (HMESH_MPI)
               CALL computeMPIFaceFlux ( f )

            end select
            end associate
         END DO           
!$omp end do
         
      END SUBROUTINE computeRiemannFluxes
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE computeElementInterfaceFlux(f)
         use FaceClass
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         real(kind=RP) :: inv_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
!         
!
!        --------------
!        Invscid fluxes:
!        --------------
!
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  flux   = inv_flux(:,i,j) )

               CALL ViscousMethod % RiemannSolver(QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  flux  = visc_flux(:,i,j) )
               
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j) ) * f % geom % scal(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         call f % ProjectFluxToElements(flux, (/1,2/))

      END SUBROUTINE computeElementInterfaceFlux

      SUBROUTINE computeMPIFaceFlux(f)
         use FaceClass
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: inv_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
!         
!
!        --------------
!        Invscid fluxes:
!        --------------
!
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  flux   = inv_flux(:,i,j) )

               CALL ViscousMethod % RiemannSolver(QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  flux  = visc_flux(:,i,j) )
               
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j) ) * f % geom % scal(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectFluxToElements(flux, (/thisSide, HMESH_NONE/))


      end subroutine ComputeMPIFaceFlux

      SUBROUTINE computeBoundaryFlux(f, time, externalStateProcedure , externalGradientsProcedure)
      USE ElementClass
      use FaceClass
      USE DGViscousDiscretization
      USE Physics
      USE BoundaryConditionFunctions
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(Face),    intent(inout) :: f
      REAL(KIND=RP)                :: time
      EXTERNAL                     :: externalStateProcedure
      EXTERNAL                     :: externalGradientsProcedure
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER, DIMENSION(2)           :: N
      REAL(KIND=RP)                   :: inv_flux(N_EQN)
      REAL(KIND=RP)                   :: UGradExt(NDIM , N_GRAD_EQN) , visc_flux(N_EQN)
      real(kind=RP)                   :: fStar(N_EQN, 0:f % Nf(1), 0: f % Nf(2))
      CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryType
            
      boundaryType = f % boundaryType
      
      DO j = 0, f % Nf(2)
         DO i = 0, f % Nf(1)
!
!           Inviscid part
!           -------------
            f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
            CALL externalStateProcedure( f % geom % x(:,i,j), &
                                         time, &
                                         f % geom % normal(:,i,j), &
                                         f % storage(2) % Q(:,i,j),&
                                         boundaryType )
            CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                               QRight = f % storage(2) % Q(:,i,j), &   
                               nHat   = f % geom % normal(:,i,j), &
                               flux   = inv_flux)
!
!           ViscousPart
!           -----------
            if ( flowIsNavierStokes ) then

            
            UGradExt(IX,:) = f % storage(1) % U_x(:,i,j)
            UGradExt(IY,:) = f % storage(1) % U_y(:,i,j)
            UGradExt(IZ,:) = f % storage(1) % U_z(:,i,j)

            CALL externalGradientsProcedure(  f % geom % x(:,i,j), &
                                              time, &
                                              f % geom % normal(:,i,j), &
                                              UGradExt,&
                                              boundaryType )

            CALL ViscousMethod % RiemannSolver( QLeft = f % storage(1) % Q(:,i,j), &
                                                QRight = f % storage(2) % Q(:,i,j), &
                                                U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                U_xRight = UGradExt(IX,:) , &
                                                U_yRight = UGradExt(IY,:) , &
                                                U_zRight = UGradExt(IZ,:) , &
                                                nHat = f % geom % normal(:,i,j), &
                                                flux = visc_flux )
            else
               visc_flux = 0.0_RP
            end if

            fStar(:,i,j) = (inv_flux - visc_flux) * f % geom % scal(i,j)
         END DO   
      END DO   

      call f % ProjectFluxToElements(fStar, (/1, HMESH_NONE/))

      END SUBROUTINE computeBoundaryFlux
!
!//////////////////////////////////////////////////////////////////////// 
! 
      REAL(KIND=RP) FUNCTION MaxTimeStep( self, cfl ) 
         IMPLICIT NONE
         TYPE(DGSem)    :: self
         REAL(KIND=RP)  :: cfl
         
         
         MaxTimeStep  = cfl/MaximumEigenvalue( self )
      
      END FUNCTION MaxTimeStep
!
!////////////////////////////////////////////////////////////////////////
!
   FUNCTION MaximumEigenvalue( self )
!
!  -------------------------------------------------------------------
!  Estimate the maximum eigenvalue of the system. This
!  routine computes a heuristic based on the smallest mesh
!  spacing (which goes as 1/N^2) AND the eigenvalues of the
!  particular hyperbolic system being solved. This is not
!  the only way to estimate the eigenvalues, but it works in practice.
!  Other people use variations on this and we make no assertions that
!  it is the only or best way. Other variations look only at the smallest
!  mesh values, other account for differences across the element.
!  -------------------------------------------------------------------
!
      USE Physics
      use MPI_Process_Info
#ifdef _HAS_MPI_
      use mpi
#endif
      
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(DGSem), TARGET    :: self
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP)               :: eValues(3)
      REAL(KIND=RP)               :: dcsi, deta, dzet, jac, lamcsi, lamzet, lameta
      INTEGER                     :: i, j, k, id
      INTEGER                     :: N(3)
      REAL(KIND=RP)               :: MaximumEigenvalue
      real(kind=RP)               :: localMaximumEigenvalue
      EXTERNAL                    :: ComputeEigenvaluesForState
      REAL(KIND=RP)               :: Q(N_EQN)
      TYPE(NodalStorage), POINTER :: spAxi_p, spAeta_p, spAzeta_p
!            
      
!
!     -----------------------------------------------------------
!     TODO:
!     dcsi, deta and dzet have been modified so they work for N=2
!     However, I'm not sure if this modification is OK.
!     Besides, problems are expected for N=0.
!     -----------------------------------------------------------
!
      MaximumEigenvalue = 0.0_RP
!$omp parallel shared(self,MaximumEigenvalue) default(private) 
!$omp do reduction(max:MaximumEigenvalue)
      DO id = 1, SIZE(self % mesh % elements) 
         N = self % mesh % elements(id) % Nxyz
         spAxi_p => self % spA(N(1))
         spAeta_p => self % spA(N(2))
         spAzeta_p => self % spA(N(3))
         IF ( ANY(N<1) ) THEN 
            PRINT*, "Error in MaximumEigenvalue function (N<1)"    
         ENDIF         
         
         dcsi = 1.0_RP / abs( spAxi_p   % x(1) - spAxi_p   % x (0) )   
         deta = 1.0_RP / abs( spAeta_p  % x(1) - spAeta_p  % x (0) )
         dzet = 1.0_RP / abs( spAzeta_p % x(1) - spAzeta_p % x (0) )
         DO k = 0, N(3)
            DO j = 0, N(2)
               DO i = 0, N(1)
!
!                 ------------------------------------------------------------
!                 The maximum eigenvalues for a particular state is determined
!                 by the physics.
!                 ------------------------------------------------------------
!
                  Q(1:N_EQN) = self % mesh % elements(id) % storage % Q(1:N_EQN,i,j,k)
                  CALL ComputeEigenvaluesForState( Q , eValues )
!
!                 ----------------------------
!                 Compute contravariant values
!                 ----------------------------
!              
                  lamcsi =  ( self % mesh % elements(id) % geom % jGradXi(IX,i,j,k)   * eValues(IX) + &
        &                     self % mesh % elements(id) % geom % jGradXi(IY,i,j,k)   * eValues(IY) + &
        &                     self % mesh % elements(id) % geom % jGradXi(IZ,i,j,k)   * eValues(IZ) ) * dcsi
        
                  lameta =  ( self % mesh % elements(id) % geom % jGradEta(IX,i,j,k)  * eValues(IX) + &
        &                     self % mesh % elements(id) % geom % jGradEta(IY,i,j,k)  * eValues(IY) + &
        &                     self % mesh % elements(id) % geom % jGradEta(IZ,i,j,k)  * eValues(IZ) ) * deta
        
                  lamzet =  ( self % mesh % elements(id) % geom % jGradZeta(IX,i,j,k) * eValues(IX) + &
        &                     self % mesh % elements(id) % geom % jGradZeta(IY,i,j,k) * eValues(IY) + &
        &                     self % mesh % elements(id) % geom % jGradZeta(IZ,i,j,k) * eValues(IZ) ) * dzet
        
                  jac               = self % mesh % elements(id) % geom % jacobian(i,j,k)
                  MaximumEigenvalue = MAX( MaximumEigenvalue, ABS(lamcsi/jac) + ABS(lameta/jac) + ABS(lamzet/jac) )
               END DO
            END DO
         END DO
      END DO 
!$omp end do
!$omp end parallel
      
   END FUNCTION MaximumEigenvalue


   END Module DGSEMClass
