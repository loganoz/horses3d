!
!//////////////////////////////////////////////////////
!
!   @File:    DGSEMClass.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Jan 14 17:14:37 2018
!   @Last revision date: Wed Jan 31 18:27:04 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 1181c365aba00e78739d327d06901d6d8ca99e02
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
Module DGSEMClass
   USE NodalStorageClass
   USE HexMeshClass
   USE PhysicsStorage
   USE SpatialDiscretization
   use MonitorsClass
   use Physics
#ifdef _HAS_MPI_
   use mpi
#endif
   
   IMPLICIT NONE

   private
   public DGSem, ComputeTimeDerivative, ComputeMaxResidual, InitializeNodalStorage
   public DestructGlobalNodalStorage
   
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
      integer                                                 :: nodes                 ! Either GAUSS or GAUSLOBATTO
      INTEGER                                                 :: numberOfTimeSteps
      INTEGER                                                 :: NDOF                         ! Number of degrees of freedom
      INTEGER           , ALLOCATABLE                         :: Nx(:), Ny(:), Nz(:)
      TYPE(HexMesh)                                           :: mesh
      PROCEDURE(externalStateSubroutine)    , NOPASS, POINTER :: externalState     => NULL()
      PROCEDURE(externalGradientsSubroutine), NOPASS, POINTER :: externalGradients => NULL()
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
      integer                     :: dir2D
      character(len=LINE_LENGTH)  :: meshFileName
      character(len=*), parameter :: TWOD_OFFSET_DIR_KEY = "2d mesh offset direction"
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
         call constructMeshFromFile( self % mesh, meshfileName, nodes, Nx, Ny, Nz, MeshInnerCurves , dir2D, success )
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
      CALL constructMeshFromFile( self % mesh, meshfileName, nodes, Nx, Ny, Nz, MeshInnerCurves , dir2D, success )
      
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
      call Initialize_SpaceAndTimeMethods(controlVariables, self % mesh)
      
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
      use ViscousMethods
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
!////////////////////////////////////////////////////////////////////////////////////////      
!
!  -----------------------------------
!  Compute maximum residual L_inf norm
!  -----------------------------------
   FUNCTION ComputeMaxResidual(self) RESULT(maxResidual)
      use MPI_Process_Info
      IMPLICIT NONE
      !----------------------------------------------
      CLASS(DGSem)  :: self
      REAL(KIND=RP) :: maxResidual(N_EQN)
      !----------------------------------------------
      INTEGER       :: id , eq, ierr
      real(kind=RP) :: cMax
      REAL(KIND=RP) :: localCMax(NCONS)
      !----------------------------------------------
      
      maxResidual = 0.0_RP

      cMax = 0.0_RP
      localCMax = 0.0_RP
!$omp parallel shared(maxResidual,cMax,self) default(private)
!$omp do reduction(max:cMax)
      DO id = 1, SIZE( self % mesh % elements )
         localCMax(1) = maxval(abs(self % mesh % elements(id) % storage % QDot(1,:,:,:)))
         
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
         use Physics, only: QuarticDWPDerivative
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
         INTEGER :: k, eID
         interface
            subroutine UserDefinedSourceTerm(mesh, time, thermodynamics_, dimensionless_, refValues_)
               USE HexMeshClass
               use PhysicsStorage
               IMPLICIT NONE
               CLASS(HexMesh)                        :: mesh
               REAL(KIND=RP)                         :: time
               type(Thermodynamics_t),    intent(in) :: thermodynamics_
               type(Dimensionless_t),     intent(in) :: dimensionless_
               type(RefValues_t),         intent(in) :: refValues_
            end subroutine UserDefinedSourceTerm
         end interface

!
!        **************************************
!        Compute chemical potential: Q stores c
!        **************************************
!
!$omp parallel shared(self, time)
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call self % mesh % ProlongSolutionToFaces()
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
!$omp single
         call self % mesh % UpdateMPIFacesSolution
!$omp end single
!
!        -----------------
!        Compute gradients
!        -----------------
!
         CALL DGSpatial_ComputeGradient( self % mesh , time , self % externalState , self % externalGradients )

!$omp single
         call self % mesh % UpdateMPIFacesGradients
!$omp end single
!
!        ------------------------------
!        Compute the chemical potential
!        ------------------------------
!
         call ComputeLaplacian(mesh = self % mesh , &
                               t    = time, &
                  externalState     = self % externalState, &
                  externalGradients = self % externalGradients )

         associate(c_alpha => thermodynamics % c_alpha, &
                   c_beta  => thermodynamics % c_beta    ) 
!$omp do
         do eID = 1, self % mesh % no_of_elements
            associate(e => self % mesh % elements(eID))
            e % storage % mu = - POW2(dimensionless % eps) * e % storage % QDot(1,:,:,:)
            e % storage % c  = e % storage % Q(1,:,:,:)
            call QuarticDWPDerivative(e % Nxyz, e % storage % c, c_alpha, c_beta, e % storage % mu)
            e % storage % Q(1,:,:,:) = e % storage % mu
            end associate
         end do
!$omp end do
         end associate
!
!        *************************
!        Compute cDot: Q stores mu
!        *************************
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call self % mesh % ProlongSolutionToFaces()
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
!$omp single
         call self % mesh % UpdateMPIFacesSolution
!$omp end single
!
!        -----------------
!        Compute gradients
!        -----------------
!
         CALL DGSpatial_ComputeGradient( self % mesh , time , self % externalState , self % externalGradients )
!
!$omp single
         call self % mesh % UpdateMPIFacesGradients
!$omp end single
!
!        ------------------------------
!        Compute the chemical potential
!        ------------------------------
!
         call ComputeLaplacian(mesh = self % mesh , &
                               t    = time, &
                  externalState     = self % externalState, &
                  externalGradients = self % externalGradients )
!
!        Add a source term
!        -----------------
!$omp single
         call UserDefinedSourceTerm(self % mesh, time, thermodynamics, dimensionless, refValues)
!$omp end single
!
!        *****************************
!        Return the concentration to Q
!        *****************************
!
!$omp do
         do eID = 1, self % mesh % no_of_elements
            associate(e => self % mesh % elements(eID))
            e % storage % Q(1,:,:,:) = e % storage % c
            end associate
         end do
!$omp end do
!$omp end parallel

      END SUBROUTINE ComputeTimeDerivative
!
!////////////////////////////////////////////////////////////////////////
!
end module DGSEMClass
