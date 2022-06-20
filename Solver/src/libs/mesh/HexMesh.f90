!
!//////////////////////////////////////////////////////
!
!   @File:
!   @Last revision commit: db0669408a09b12909b03856ec101a38d42f4e79
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
MODULE HexMeshClass
      use Utilities                       , only: toLower, almostEqual, AlmostEqualRelax
      use SMConstants
      USE MeshTypes
      USE NodeClass
      USE ElementClass
      USE FaceClass
      use FacePatchClass
      USE TransfiniteMapClass
      use SharedBCModule
      use ElementConnectivityDefinitions
      use ZoneClass                       , only: Zone_t, ConstructZones, ReassignZones
      use PhysicsStorage
      use NodalStorageClass
      use MPI_Process_Info
      use MPI_Face_Class
      use FluidData
      use StorageClass
      use FileReadingUtilities            , only: RemovePath, getFileName
      use FTValueDictionaryClass          , only: FTValueDictionary
      use SolutionFile
      use BoundaryConditions,               only: BCs
      use IntegerDataLinkedList           , only: IntegerDataLinkedList_t
      use PartitionedMeshClass            , only: mpi_partition
      use IBMClass
#if defined(NAVIERSTOKES)
      use WallDistance
#endif
#ifdef _HAS_MPI_
      use mpi
#endif
      IMPLICIT NONE

      private
      public      HexMesh
      public      Neighbor_t, NUM_OF_NEIGHBORS

      public      GetOriginalNumberOfFaces
      public      ConstructFaces, ConstructPeriodicFaces
      public      DeletePeriodicMinusFaces, GetElementsFaceIDs
      public      no_of_stats_variables
!
!     ---------------
!     Mesh definition
!     ---------------
!
      type HexMesh
         integer                                   :: numberOfFaces
         integer                                   :: nodeType
         integer                                   :: no_of_elements
         integer                                   :: no_of_allElements
         integer                                   :: no_of_faces
         integer                                   :: dt_restriction = DT_FIXED     ! Time step restriction of last step (DT_FIXED -initial value-, DT_DIFF or DT_CONV)
         integer      , dimension(:), allocatable  :: Nx, Ny, Nz
         integer                                   :: NDOF
         integer,                     allocatable  :: faces_interior(:)
         integer,                     allocatable  :: faces_mpi(:)
         integer,                     allocatable  :: faces_boundary(:)
         integer,                     allocatable  :: elements_sequential(:)
         integer,                     allocatable  :: elements_mpi(:)
         integer, allocatable                      :: HOPRnodeIDs(:)
         character(len=LINE_LENGTH)                :: meshFileName
         type(SolutionStorage_t)                   :: storage              ! Here the solution and its derivative are stored
         type(Node)   , dimension(:), allocatable  :: nodes
         type(Face)   , dimension(:), allocatable  :: faces
         type(Element), dimension(:), allocatable  :: elements
         type(MPI_FacesSet_t)                      :: MPIfaces
         type(IBM_type)                            :: IBM
         class(Zone_t), dimension(:), allocatable  :: zones
         logical                                   :: child       = .FALSE.         ! Is this a (multigrid) child mesh? default .FALSE.
         logical                                   :: meshIs2D    = .FALSE.         ! Is this a 2D mesh? default .FALSE.
         integer                                   :: dir2D                         ! If it is in fact a 2D mesh, dir 2D stores the global direction IX, IY or IZ
         integer                                   :: dir2D_ctrl                    ! dir2D as in the control file
         logical                                   :: anisotropic = .FALSE.         ! Is the mesh composed by elements with anisotropic polynomial orders? default false
         logical                                   :: ignoreBCnonConformities = .FALSE.
         contains
            procedure :: destruct                      => HexMesh_Destruct
            procedure :: Describe                      => HexMesh_Describe
            procedure :: DescribePartition             => DescribeMeshPartition
            procedure :: AllocateStorage               => HexMesh_AllocateStorage
            procedure :: ConstructZones                => HexMesh_ConstructZones
            procedure :: DefineAsBoundaryFaces         => HexMesh_DefineAsBoundaryFaces
            procedure :: CheckIfMeshIs2D               => HexMesh_CheckIfMeshIs2D
            procedure :: CorrectOrderFor2DMesh         => HexMesh_CorrectOrderFor2DMesh
            procedure :: SetConnectivitiesAndLinkFaces => HexMesh_SetConnectivitiesAndLinkFaces
            procedure :: UpdateFacesWithPartition      => HexMesh_UpdateFacesWithPartition
            procedure :: ConstructGeometry             => HexMesh_ConstructGeometry
            procedure :: ProlongSolutionToFaces        => HexMesh_ProlongSolutionToFaces
            procedure :: ProlongGradientsToFaces       => HexMesh_ProlongGradientsToFaces
            procedure :: PrepareForIO                  => HexMesh_PrepareForIO
            procedure :: Export                        => HexMesh_Export
            procedure :: ExportOrders                  => HexMesh_ExportOrders
            procedure :: ExportBoundaryMesh            => HexMesh_ExportBoundaryMesh
            procedure :: SaveSolution                  => HexMesh_SaveSolution
            procedure :: pAdapt                        => HexMesh_pAdapt
#if defined(NAVIERSTOKES)
            procedure :: SaveStatistics                => HexMesh_SaveStatistics
            procedure :: ResetStatistics               => HexMesh_ResetStatistics
#endif
            procedure :: LoadSolution                  => HexMesh_LoadSolution
            procedure :: LoadSolutionForRestart        => HexMesh_LoadSolutionForRestart
            procedure :: WriteCoordFile
            procedure :: UpdateMPIFacesPolynomial      => HexMesh_UpdateMPIFacesPolynomial
            procedure :: UpdateMPIFacesSolution        => HexMesh_UpdateMPIFacesSolution
            procedure :: UpdateMPIFacesGradients       => HexMesh_UpdateMPIFacesGradients
            procedure :: UpdateMPIFacesAviscflux       => HexMesh_UpdateMPIFacesAviscflux
            procedure :: GatherMPIFacesSolution        => HexMesh_GatherMPIFacesSolution
            procedure :: GatherMPIFacesGradients       => HexMesh_GatherMPIFacesGradients
            procedure :: GatherMPIFacesAviscFlux       => HexMesh_GatherMPIFacesAviscFlux
            procedure :: FindPointWithCoords           => HexMesh_FindPointWithCoords
            procedure :: FindPointWithCoordsInNeighbors=> HexMesh_FindPointWithCoordsInNeighbors
            procedure :: ComputeWallDistances          => HexMesh_ComputeWallDistances
            procedure :: ConformingOnZone              => HexMesh_ConformingOnZone
            procedure :: SetStorageToEqn          => HexMesh_SetStorageToEqn
#if defined(INCNS) && defined(CAHNHILLIARD)
            procedure :: ConvertDensityToPhaseFIeld    => HexMesh_ConvertDensityToPhaseField
            procedure :: ConvertPhaseFieldToDensity    => HexMesh_ConvertPhaseFieldToDensity
#endif
            procedure :: copy                          => HexMesh_Assign
            generic   :: assignment(=)                 => copy
      end type HexMesh

      integer, parameter :: NUM_OF_NEIGHBORS = 6 ! Hardcoded: Hexahedral conforming meshes
      integer            :: no_of_stats_variables
      
      TYPE Neighbor_t         ! added to introduce colored computation of numerical Jacobian (is this the best place to define this type??) - only usable for conforming meshes
         INTEGER :: elmnt(NUM_OF_NEIGHBORS+1) ! "7" hardcoded for 3D hexahedrals in conforming meshes (the last one is itself)... This definition must change if the code is expected to be more general
      END TYPE Neighbor_t
!
!     ========
      CONTAINS
!     ========
!
!     -----------
!     Destructors
!     -----------
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE HexMesh_Destruct( self )
         IMPLICIT NONE
         CLASS(HexMesh) :: self

         safedeallocate (self % Nx)
         safedeallocate (self % Ny)
         safedeallocate (self % Nz)
!
!        -----
!        Nodes
!        -----
!
         call self % nodes % destruct
         DEALLOCATE( self % nodes )
         safedeallocate (self % HOPRnodeIDs)
!
!        --------
!        Elements
!        --------
!
         call self % elements % destruct
         DEALLOCATE( self % elements )
!
!        -----
!        Faces
!        -----
!
         call self % faces % Destruct
         DEALLOCATE( self % faces )
!
!        -----
!        Zones
!        -----
!
         if (allocated(self % zones)) DEALLOCATE( self % zones )
!
!        ----------------
!        Solution storage
!        ----------------
!
         call self % storage % destruct

         safedeallocate(self % elements_sequential)
         safedeallocate(self % elements_mpi)
         safedeallocate(self % faces_interior)
         safedeallocate(self % faces_mpi)
         safedeallocate(self % faces_boundary)

!
!        ----------------
!        IBM storage
!        ----------------
!
         if( self% IBM% active ) then
            if( self% child ) then
               call self% IBM% destruct( .true. )
            else
               call self% IBM% destruct( .false. )
            end if
         end if
         
      END SUBROUTINE HexMesh_Destruct
!
!     -------------
!     Print methods
!     -------------
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PrintMesh( self )
      IMPLICIT NONE
      TYPE(HexMesh) :: self
      INTEGER ::  k

      PRINT *, "Nodes..."
      DO k = 1, SIZE(self % nodes)
         CALL PrintNode( self % nodes(k), k )
      END DO
      PRINT *, "Elements..."
      DO k = 1, SIZE(self % elements)
         CALL PrintElement( self % elements(k), k )
      END DO
      PRINT *, "Faces..."
      DO k = 1, SIZE(self % faces)
         CALL  self % faces(k) % Print
      END DO

      END SUBROUTINE PrintMesh
!
!////////////////////////////////////////////////////////////////////////
!
      integer function GetOriginalNumberOfFaces(self)
         USE FTMultiIndexTableClass
         USE FTValueClass

         IMPLICIT NONE
         TYPE(HexMesh) :: self

         INTEGER                 :: eID, faceNumber
         INTEGER                 :: faceID
         INTEGER                 :: nodeIDs(8), faceNodeIDs(4), j

         CLASS(FTMultiIndexTable), POINTER  :: table
         CLASS(FTObject), POINTER :: obj
         CLASS(FTValue) , POINTER :: v

         ALLOCATE(table)
         CALL table % initWithSize( SIZE( self % nodes) )

         GetOriginalNumberOfFaces = 0
         DO eID = 1, SIZE( self % elements )

            nodeIDs = self % elements(eID) % nodeIDs
            DO faceNumber = 1, 6
               DO j = 1, 4
                  faceNodeIDs(j) = nodeIDs(localFaceNode(j,faceNumber))
               END DO

               IF (.not. table % containsKeys(faceNodeIDs) )     THEN
                  GetOriginalNumberOfFaces = GetOriginalNumberOfFaces + 1

                  ALLOCATE(v)
                  CALL v % initWithValue(GetOriginalNumberOfFaces)
                  obj => v
                  CALL table % addObjectForKeys(obj,faceNodeIDs)
                  CALL release(v)
               END IF
            END DO
         END DO

         CALL release(table)

      end function GetOriginalNumberOfFaces

      SUBROUTINE ConstructFaces( self, success )
!
!     -------------------------------------------------------------
!     Go through the elements and find the unique faces in the mesh
!     -------------------------------------------------------------
!
         use IntegerArrayLinkedListTable
         IMPLICIT NONE
         TYPE(HexMesh) :: self
         LOGICAL       :: success

         INTEGER                 :: eID, faceNumber
         INTEGER                 :: faceID
         INTEGER                 :: nodeIDs(8), faceNodeIDs(4), j
         type(Table_t)           :: table

         table = Table_t(size(self % nodes))

         self % numberOfFaces = 0
         DO eID = 1, SIZE( self % elements )

            nodeIDs = self % elements(eID) % nodeIDs
            DO faceNumber = 1, 6
               DO j = 1, 4
                  faceNodeIDs(j) = nodeIDs(localFaceNode(j,faceNumber))
               END DO

               faceID = table % ContainsEntry(faceNodeIDs)
               IF ( faceID .ne. 0 )     THEN
!
!                 --------------------------------------------------------------
!                 Add this element to the slave side of the face associated with
!                 these nodes.
!                 --------------------------------------------------------------
!
                  self % faces(faceID) % elementIDs(2)  = eID
                  self % faces(faceID) % elementSide(2) = faceNumber
                  self % faces(faceID) % FaceType       = HMESH_INTERIOR
                  self % faces(faceID) % rotation       = faceRotation(masterNodeIDs = self % faces(faceID) % nodeIDs, &
                                                                       slaveNodeIDs  = faceNodeIDs                      )
               ELSE
!
!                 ------------------
!                 Construct new face
!                 ------------------
!
                  self % numberOfFaces = self % numberOfFaces + 1

                  IF(self % numberOfFaces > SIZE(self % faces))     THEN

                     call table % Destruct
                     PRINT *, "Too many faces for # of elements:", self % numberOfFaces, " vs ", SIZE(self % faces)
                     success = .FALSE.
                     RETURN
                  END IF

                  CALL self % faces(self % numberOfFaces) % Construct(ID  = self % numberOfFaces, &
                                                                      nodeIDs = faceNodeIDs, &
                                                                      elementID = eID,       &
                                                                      side = faceNumber)

                  self % faces(self % numberOfFaces) % boundaryName = &
                           self % elements(eID) % boundaryName(faceNumber)
!
!                 ----------------------------------------------
!                 Mark which face is associated with these nodes
!                 ----------------------------------------------
!
                  call table % AddEntry(faceNodeIDs)
               END IF
            END DO

         END DO

         call table % Destruct

      END SUBROUTINE ConstructFaces

      subroutine GetElementsFaceIDs(self)
         implicit none
         type(HexMesh), intent(inout)  :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: fID, eL, eR, e, side

         do fID = 1, size(self % faces)
            select case (self % faces(fID) % faceType)
            case (HMESH_INTERIOR)
               eL = self % faces(fID) % elementIDs(1)
               eR = self % faces(fID) % elementIDs(2)

               self % elements(eL) % faceIDs(self % faces(fID) % elementSide(1)) = fID
               self % elements(eR) % faceIDs(self % faces(fID) % elementSide(2)) = fID
               self % elements(eL) % faceSide(self % faces(fID) % elementSide(1)) = 1
               self % elements(eR) % faceSide(self % faces(fID) % elementSide(2)) = 2

            case (HMESH_BOUNDARY)
               eL = self % faces(fID) % elementIDs(1)
               self % elements(eL) % faceIDs(self % faces(fID) % elementSide(1)) = fID
               self % elements(eL) % faceSide(self % faces(fID) % elementSide(1)) = 1

            case (HMESH_MPI)
               side = maxloc(self % faces(fID) % elementIDs, 1)

               e = self % faces(fID) % elementIDs(side)
               self % elements(e) % faceIDs(self % faces(fID) % elementSide(side)) = fID
               self % elements(e) % faceSide(self % faces(fID) % elementSide(side)) = side
               self % elements(e) % hasSharedFaces = .true.

            case (HMESH_UNDEFINED)
               eL = self % faces(fID) % elementIDs(1)
               self % elements(eL) % faceIDs(self % faces(fID) % elementSide(1)) = fID
               self % elements(eL) % faceSide(self % faces(fID) % elementSide(1)) = 1

            end select
         end do

      end subroutine GetElementsFaceIDs
!
!////////////////////////////////////////////////////////////////////////
!
!
!---------------------------------------------------------------------
!! Element faces can be rotated with respect to each other. Orientation
!! gives the relative orientation of the master (1) face to the
!! slave (2) face . In this routine,
!! orientation is measured in 90 degree increments:
!!                   rotation angle = orientation*pi/2
!!
!! As an example, faceRotation = 1 <=> rotating master by 90 deg.
!
      INTEGER pure FUNCTION faceRotation(masterNodeIDs, slaveNodeIDs)
         IMPLICIT NONE
         INTEGER, DIMENSION(4), intent(in) :: masterNodeIDs, slaveNodeIDs !< Node IDs
!
!        ---------------
!        Local variables
!        ---------------
!
         integer, dimension(4), parameter :: NEXTNODE = (/2,3,4,1/)
         INTEGER :: j
!
!        Rotate until both first nodes match (each j corresponds to a 90deg rotation)
!        -----------------------------------
         DO j = 1, 4
            IF(masterNodeIDs(1) == slaveNodeIDs(j)) EXIT
         END DO
!
!        Check whether the orientation is same or opposite
!        -------------------------------------------------
         if ( masterNodeIDS(2) == slaveNodeIDs(NEXTNODE(j)) ) then
            faceRotation = j - 1
         else
            faceRotation = j + 3
         end if

      END FUNCTION faceRotation
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructPeriodicFaces(self, useRelaxTol)
      USE Physics
      IMPLICIT NONE
!
!-------------------------------------------------------------------
! This subroutine looks for periodic boundary conditions. If they
! are found, periodic+ face is set as an interior face. The slave
! face is the periodic- face and will be deleted in the following
! step.
!-------------------------------------------------------------------
!
!
!--------------------
! External variables
!--------------------
!
      TYPE(HexMesh)              :: self
      LOGICAL, intent(in)        :: useRelaxTol

!
!--------------------
! Local variables
!--------------------
!
!
      REAL(KIND=RP)              :: x1(NDIM), x2(NDIM), edge_length(4), min_edge_length
      LOGICAL                    :: master_matched(4), slave_matched(4), success, found
      INTEGER                    :: coord, slaveNodeIDs(4), localCoord

      INTEGER                    :: i,j,k,l
      integer                    :: zIDplus, zIDMinus, iFace, jFace
      character(len=LINE_LENGTH) :: associatedBname
!
!     --------------------------------------------
!     Loop to find faces with the label "periodic"
!     --------------------------------------------
!
!     -----------------------------
!     Loop zones with BC "periodic"
!     -----------------------------
!

      do zIDPlus = 1, size(self % zones)
!
!        Cycle if the zone is not periodic
!        ---------------------------------
         if ( trim(BCs(zIDPlus) % bc % bcType) .ne. "periodic") cycle
!
!        Cycle if the zone has already marked to be deleted
!        --------------------------------------------------
         if ( self % zones(zIDPlus) % toBeDeleted ) cycle
!
!        Reset the coordinate (changes when changing zones)
!        --------------------------------------------------
         coord = 0
!
!        Get the marker of the associated zone
!        -------------------------------------
         found = .false.
         do zIDMinus = 1, size(self % zones)
            call BCs(zIDPlus) % bc % GetPeriodicPair(associatedBname)
            if ( trim(associatedBname) .eq. trim(self % zones(zIDMinus) % Name) ) then
               found = .true.
               self % zones(zIDMinus) % toBeDeleted = .true.
               exit
            end if
         end do

         if ( .not. found ) then
            print*, 'coupled boundary "',trim(associatedBname),' for boundary "',trim(self % zones(zIDPlus) % Name),'" not found.'
            errorMessage(STD_OUT)
            stop
         end if
!
!        Loop faces in the periodic+ zone
!        --------------------------------
ploop:   do iFace = 1, self % zones(zIDPlus) % no_of_faces
!
!           Get the face ID
!           ---------------
            i = self % zones(zIDPlus) % faces(iFace)
!
!           Consider only HMESH_UNDEFINED faces
!           -----------------------------------
            if ( (self % faces(i) % faceType .ne. HMESH_UNDEFINED)) cycle ploop
!
!           Loop faces in the periodic- zone
!           --------------------------------
mloop:      do jFace = 1, self % zones(zIDMinus) % no_of_faces
!
!              Get the face ID
!              ---------------
               j = self % zones(zIDMinus) % faces(jFace)
!
!              Consider only HMESH_UNDEFINED faces
!              -----------------------------------
               if ( (self % faces(j) % faceType .ne. HMESH_UNDEFINED)) cycle mloop
!
!              ----------------------------------------------------------------------------------------
!              The index i is a periodic+ face
!              The index j is a periodic- face
!              We are looking for couples of periodic+ and periodic- faces where 2 of the 3 coordinates
!              in all the corners are shared. The non-shared coordinate has to be always the same one.
!              ---------------------------------------------------------------------------------------
!
               master_matched(:)   = .FALSE.     ! True if the master corner finds a partner
               slave_matched(:)    = .FALSE.     ! True if the slave corner finds a partner

               ! compute minimum edge length to make matching tolerance relative to element size. Assumes that the periodic
               ! coordinate of all the nodes not vary significativelly compared to the other two coordinates.
               if (useRelaxTol) then
                   edge_length(1)=NORM2(self % nodes(self % faces(i) % nodeIDs(1)) % x - self % nodes(self % faces(i) % nodeIDs(2)) % x)
                   edge_length(2)=NORM2(self % nodes(self % faces(i) % nodeIDs(2)) % x - self % nodes(self % faces(i) % nodeIDs(3)) % x)
                   edge_length(3)=NORM2(self % nodes(self % faces(i) % nodeIDs(3)) % x - self % nodes(self % faces(i) % nodeIDs(4)) % x)
                   edge_length(4)=NORM2(self % nodes(self % faces(i) % nodeIDs(4)) % x - self % nodes(self % faces(i) % nodeIDs(1)) % x)
                   min_edge_length=minval(edge_length)
               end if

               if ( coord .eq. 0 ) then
!
!                 Check all coordinates
!                 ---------------------
                  do localCoord = 1, 3
                     master_matched = .false.
                     slave_matched = .false.
mastercoord:         DO k = 1, 4
                        x1 = self%nodes(self%faces(i)%nodeIDs(k))%x
slavecoord:             DO l = 1, 4
                           IF (.NOT.slave_matched(l)) THEN
                              x2 = self%nodes(self%faces(j)%nodeIDs(l))%x
                              IF (useRelaxTol) THEN
                                  CALL CompareTwoNodesRelax(x1, x2, master_matched(k), localCoord, min_edge_length)
                              ELSE
                                  CALL CompareTwoNodes(x1, x2, master_matched(k), localCoord)
                              END IF
                              IF (master_matched(k)) THEN
                                 slave_matched(l) = .TRUE.
                                 EXIT  slavecoord
                              ENDIF
                           ENDIF
                        ENDDO    slavecoord
                        IF (.NOT.master_matched(k)) EXIT mastercoord
                     ENDDO mastercoord

                     if ( all(master_matched) ) exit
                  end do

               else
!
!                 Check only the shared coordinates
!                 ---------------------------------
                  DO k = 1, 4
                     x1 = self%nodes(self%faces(i)%nodeIDs(k))%x
                     DO l = 1, 4
                        IF (.NOT.slave_matched(l)) THEN
                           x2 = self%nodes(self%faces(j)%nodeIDs(l))%x
                           IF (useRelaxTol) THEN
                               CALL CompareTwoNodesRelax(x1, x2, master_matched(k), coord, min_edge_length)
                           ELSE
                               CALL CompareTwoNodes(x1, x2, master_matched(k), coord)
                           END IF
                           IF (master_matched(k)) THEN
                              slave_matched(l) = .TRUE.
                              EXIT
                           ENDIF
                        ENDIF
                     ENDDO
                     IF (.NOT.master_matched(k)) EXIT
                  ENDDO

               end if

               IF ( all(master_matched) ) THEN
                  if ( coord .eq. 0 ) coord = localCoord
                  self % faces(i) % boundaryName   = emptyBCName
                  self % faces(i) % elementIDs(2)  = self % faces(j) % elementIDs(1)
                  self % faces(i) % elementSide(2) = self % faces(j) % elementSide(1)
                  self % faces(i) % FaceType       = HMESH_INTERIOR
                  self % elements(self % faces(i) % elementIDs(1)) % boundaryName(self % faces(i) % elementSide(1)) = emptyBCName
                  self % elements(self % faces(i) % elementIDs(2)) % boundaryName(self % faces(i) % elementSide(2)) = emptyBCName
!
!                 To obtain the face rotation, we traduce the right element node IDs to the left
!                 ------------------------------------------------------------------------------
                  do k = 1, 4
                     x1 = self % nodes ( self % faces(i) % nodeIDs(k)) % x
                     do l = 1, 4
                        x2 = self % nodes ( self % faces(j) % nodeIDs(l) ) % x
                        IF (useRelaxTol) THEN
                            CALL CompareTwoNodesRelax(x1, x2, success, coord, min_edge_length)
                        ELSE
                            CALL CompareTwoNodes(x1, x2, success, coord)
                        END IF
                        if ( success ) then
                           slaveNodeIDs(l) = self % faces(i) % nodeIDs(k)
                        end if
                     end do
                  end do
                  self % faces(i) % rotation = faceRotation(self % faces(i) % nodeIDs, &
                                                            slaveNodeIDs)
                  cycle ploop

               ENDIF
            end do   mloop ! periodic- faces
!
!           If the code arrives here, the periodic+ face was not able to find a partner
!           ---------------------------------------------------------------------------
            print*, "When constructing periodic boundary conditions,"
            write(STD_OUT,'(A,I0,A,I0,A,I0)') "Face ",i," in zone ",zIDPlus, &
                  " was not able to find a partner. Element: ", self % faces(i) % elementIDs(1)
            errorMessage(STD_OUT)
            stop

            end do   ploop    ! periodic+ faces
         end do               ! periodic+ zones

         if ( MPI_Process % isRoot .and. useRelaxTol) print *, "Success: when matching all periodic boundary conditions with relaxed comparison"

      END SUBROUTINE ConstructPeriodicFaces
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE CompareTwoNodes(x1, x2, success, coord)
      IMPLICIT NONE
!
!-------------------------------------------------------------------
! Comparison of two nodes. If two of the three coordinates are the
! same, there is success. If there is success, the coordinate which
! is not the same is saved. If the initial value of coord is not 0,
! only that coordinate is checked.
!-------------------------------------------------------------------
!
!
!     --------------------
!     External variables
!     --------------------
!
      REAL(KIND=RP) :: x1(3)
      REAL(KIND=RP) :: x2(3)
      LOGICAL       :: success
      INTEGER       :: coord
!
!     --------------------
!     Local variables
!     --------------------
!
      INTEGER :: i
      INTEGER :: counter

      counter = 0

      IF (coord == 0) THEN

         DO i = 1,3
            IF ( AlmostEqual( x1(i), x2(i) ) ) THEN
               counter = counter + 1
            ELSE
               coord = i
            ENDIF
         ENDDO

         IF (counter.ge.2) THEN
            success = .TRUE.
         ELSE
            success = .FALSE.
         ENDIF

      ELSE

         DO i = 1,3
            IF (i /= coord) THEN
               IF ( AlmostEqual( x1(i), x2(i) ) ) THEN
                  counter = counter + 1
               ENDIF
            ENDIF
         ENDDO

         IF (counter.ge.2) THEN
            success = .TRUE.
         ELSE
            success = .FALSE.
         ENDIF

      ENDIF


      END SUBROUTINE CompareTwoNodes
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE CompareTwoNodesRelax(x1, x2, success, coord, min_edge_length)
      IMPLICIT NONE
!
!-------------------------------------------------------------------
! Similar to CompareTwoNodes, but the comparison of the two nodes
! is done relaxed by the minimum edge length
!-------------------------------------------------------------------
!     --------------------
!     External variables
!     --------------------
!
      REAL(KIND=RP) :: x1(3)
      REAL(KIND=RP) :: x2(3)
      REAL(KIND=RP) :: min_edge_length
      LOGICAL       :: success
      INTEGER       :: coord
!
!     --------------------
!     Local variables
!     --------------------
!
      INTEGER :: i
      INTEGER :: counter

      counter = 0

      IF (coord == 0) THEN

         DO i = 1,3
            IF ( AlmostEqualRelax( x1(i), x2(i) , min_edge_length ) ) THEN
               counter = counter + 1
            ELSE
               coord = i
            ENDIF
         ENDDO

         IF (counter.ge.2) THEN
            success = .TRUE.
         ELSE
            success = .FALSE.
         ENDIF

      ELSE

         DO i = 1,3
            IF (i /= coord) THEN
               IF ( AlmostEqualRelax( x1(i), x2(i) , min_edge_length ) ) THEN
                  counter = counter + 1
               ENDIF
            ENDIF
         ENDDO

         IF (counter.ge.2) THEN
            success = .TRUE.
         ELSE
            success = .FALSE.
         ENDIF

      ENDIF

      END SUBROUTINE CompareTwoNodesRelax
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DeletePeriodicminusfaces(self)
      use MPI_Face_Class
      IMPLICIT NONE
!
!-------------------------------------------------------------------
! This subroutine looks for periodic boundary conditions. If they
! are found, periodic+ face is set as an interior face. The slave
! face is the periodic- face and will be deleted in the following
! step.
!-------------------------------------------------------------------
!
!
!     --------------------
!     External variables
!     --------------------
!
      TYPE(HexMesh) :: self
!
!     --------------------
!     Local variables
!     --------------------
!
      TYPE(Face),ALLOCATABLE  :: dummy_faces(:)
      INTEGER                 :: i, domain
      INTEGER                 :: iFace, numberOfFaces
      character(len=LINE_LENGTH)    :: bName
      integer                 :: newFaceID(self % numberOfFaces)
!
!     This first loop marks which faces will not be deleted
!     -----------------------------------------------------
      newFaceID = -1
      iFace = 0
      ALLOCATE( dummy_faces(self % numberOfFaces) )
      DO i = 1, self%numberOfFaces
         if ( self % faces(i) % faceType .ne. HMESH_UNDEFINED ) then
            iFace = iFace + 1
            dummy_faces(iFace) = self%faces(i)
            dummy_faces(iFace) % ID = iFace
            newFaceID(i) = iFace
         elseif (.not. self % zones(self % faces(i) % zone) % toBeDeleted) then
            iFace = iFace + 1
            dummy_faces(iFace) = self%faces(i)
            dummy_faces(iFace) % ID = iFace
            newFaceID(i) = iFace
         ENDIF
      ENDDO

      numberOfFaces = iFace

      DEALLOCATE(self%faces)
      ALLOCATE(self%faces(numberOfFaces))

      self%numberOfFaces = numberOfFaces

      DO i = 1, self%numberOfFaces
         self%faces(i) = dummy_faces(i)
      ENDDO
!
!     Update MPI face IDs
!     -------------------
      if ( (MPI_Process % doMPIAction) .and. self % MPIfaces % Constructed ) then
         do domain = 1, MPI_Process % nProcs
            do iFace = 1, self % MPIfaces % faces(domain) % no_of_faces
               self % MPIfaces % faces(domain) % faceIDs(iFace) = newFaceID(self % MPIfaces % faces(domain) % faceIDs(iFace))
            end do
         end do
      end if
!
!     Reassign zones
!     -----------------
      CALL ReassignZones(self % faces, self % zones)

      END SUBROUTINE DeletePeriodicminusfaces
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_ProlongSolutionToFaces(self, nEqn)
         implicit none
         class(HexMesh),   intent(inout)  :: self
         integer,          intent(in)     :: nEqn
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: fIDs(6)
         integer  :: eID

!$omp do schedule(runtime)
         do eID = 1, size(self % elements)
            fIDs = self % elements(eID) % faceIDs
            call self % elements(eID) % ProlongSolutionToFaces(nEqn, &
                                                               self % faces(fIDs(1)),&
                                                               self % faces(fIDs(2)),&
                                                               self % faces(fIDs(3)),&
                                                               self % faces(fIDs(4)),&
                                                               self % faces(fIDs(5)),&
                                                               self % faces(fIDs(6)) )
         end do
!$omp end do

      end subroutine HexMesh_ProlongSolutionToFaces
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_ProlongGradientsToFaces(self, nGradEqn)
         implicit none
         class(HexMesh),   intent(inout)  :: self
         integer,          intent(in)     :: nGradEqn
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: fIDs(6)
         integer  :: eID

!$omp do schedule(runtime)
         do eID = 1, size(self % elements)
            fIDs = self % elements(eID) % faceIDs
            call self % elements(eID) % ProlongGradientsToFaces(nGradEqn, &
                                                                self % faces(fIDs(1)),&
                                                                self % faces(fIDs(2)),&
                                                                self % faces(fIDs(3)),&
                                                                self % faces(fIDs(4)),&
                                                                self % faces(fIDs(5)),&
                                                                self % faces(fIDs(6)) )
         end do
!$omp end do

      end subroutine HexMesh_ProlongGradientsToFaces
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ------------------------------------------------------------------------------------
!     HexMesh_UpdateMPIFacesPolynomial:
!        Send the face polynomial orders to the duplicate faces in the neighbor partitions
!        The information sent per face is:
!              [fNxi, fNeta, eNxi, eNeta, eNzeta, eGlobID],
!        where:
!           fNxi:    Polynomial order 1 for neighbor face
!           fNeta:   Polynomial order 2 for neighbor face
!           eNxi:    Polynomial order 1 for neighbor element
!           eNeta:   Polynomial order 2 for neighbor element
!           eNzeta:  Polynomial order 3 for neighbor element
!           eGlobID: Global ID of neighbor element
!     ------------------------------------------------------------------------------------
      subroutine HexMesh_UpdateMPIFacesPolynomial(self)
         use MPI_Face_Class
         implicit none
         !-arguments----------------------------------------------------------
         class(HexMesh)         :: self
#ifdef _HAS_MPI_
         !-local-variables----------------------------------------------------
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter
         !--------------------------------------------------------------------

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************************
!        Perform the receive request
!        ***************************
!
         do domain = 1, MPI_Process % nProcs
            call self % MPIfaces % faces(domain) % RecvN(domain)
         end do
!
!        *************
!        Send solution
!        *************
!
         do domain = 1, MPI_Process % nProcs
!
!           ---------------
!           Gather solution
!           ---------------
!
            counter = 1
            if ( self % MPIfaces % faces(domain) % no_of_faces .eq. 0 ) cycle

            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate( f => self % faces(fID))
               associate( e => self % elements(maxval(f % elementIDs)) )


               self % MPIfaces % faces(domain) % Nsend(counter:counter+1  ) = e % Nxyz(axisMap(:,f % elementSide(thisSide)))
               self % MPIfaces % faces(domain) % Nsend(counter+2:counter+4) = e % Nxyz
               self % MPIfaces % faces(domain) % Nsend(counter+5)           = e % globID

               counter = counter + 6

               end associate
               end associate
            end do
!
!           -------------
!           Send solution
!           -------------
!
            call self % MPIfaces % faces(domain) % SendN(domain)
         end do
#endif
      end subroutine HexMesh_UpdateMPIFacesPolynomial
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_UpdateMPIFacesSolution(self, nEqn)
         use MPI_Face_Class
         implicit none
         class(HexMesh)         :: self
         integer,    intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter
         integer, parameter :: otherSide(2) = (/2,1/)

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************************
!        Perform the receive request
!        ***************************
!
         do domain = 1, MPI_Process % nProcs
            call self % MPIfaces % faces(domain) % RecvQ(domain, nEqn)
         end do
!
!        *************
!        Send solution
!        *************
!
         do domain = 1, MPI_Process % nProcs
!
!           ---------------
!           Gather solution
!           ---------------
!
            counter = 1
            if ( self % MPIfaces % faces(domain) % no_of_faces .eq. 0 ) cycle

            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  self % MPIfaces % faces(domain) % Qsend(counter:counter+nEqn-1) = f % storage(thisSide) % Q(:,i,j)
                  counter = counter + nEqn
               end do               ; end do
               end associate
            end do
!
!           -------------
!           Send solution
!           -------------
!
            call self % MPIfaces % faces(domain) % SendQ(domain, nEqn)
         end do
#endif
      end subroutine HexMesh_UpdateMPIFacesSolution

      subroutine HexMesh_UpdateMPIFacesGradients(self, nEqn)
         use MPI_Face_Class
         implicit none
         class(HexMesh)      :: self
         integer, intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter
         integer, parameter :: otherSide(2) = (/2,1/)

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************************
!        Perform the receive request
!        ***************************
!
         do domain = 1, MPI_Process % nProcs
            call self % MPIfaces % faces(domain) % RecvU_xyz(domain, nEqn)
         end do
!
!        ***************
!        Gather gradients
!        ***************
!
         do domain = 1, MPI_Process % nProcs
            if ( self % MPIfaces % faces(domain) % no_of_faces .eq. 0 ) cycle

            counter = 1

            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  self % MPIfaces % faces(domain) % U_xyzsend(counter:counter+nEqn-1) = f % storage(thisSide) % U_x(:,i,j)
                  counter = counter + nEqn
               end do               ; end do

               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  self % MPIfaces % faces(domain) % U_xyzsend(counter:counter+nEqn-1) = f % storage(thisSide) % U_y(:,i,j)
                  counter = counter + nEqn
               end do               ; end do

               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  self % MPIfaces % faces(domain) % U_xyzsend(counter:counter+nEqn-1) = f % storage(thisSide) % U_z(:,i,j)
                  counter = counter + nEqn
               end do               ; end do
               end associate
            end do

            call self % MPIfaces % faces(domain) % SendU_xyz(domain, nEqn)
         end do
#endif
      end subroutine HexMesh_UpdateMPIFacesGradients
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_UpdateMPIFacesAviscFlux(self, nEqn)
         use MPI_Face_Class
         implicit none
         class(HexMesh)         :: self
         integer,    intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter
         integer, parameter :: otherSide(2) = (/2,1/)

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************************
!        Perform the receive request
!        ***************************
!
         do domain = 1, MPI_Process % nProcs
            call self % MPIfaces % faces(domain) % RecvAviscFlux(domain, nEqn)
         end do
!
!        ***********
!        Send H flux
!        ***********
!
         do domain = 1, MPI_Process % nProcs
!
!           ---------------
!           Gather solution
!           ---------------
!
            counter = 1
            if ( self % MPIfaces % faces(domain) % no_of_faces .eq. 0 ) cycle

            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  self % MPIfaces % faces(domain) % AviscFluxSend(counter:counter+nEqn-1) = &
                     f % storage(thisSide) % AviscFlux(:,i,j)
                  counter = counter + nEqn
               end do               ; end do
               end associate
            end do
!
!           -------------
!           Send solution
!           -------------
!
            call self % MPIfaces % faces(domain) % SendAviscFlux(domain, nEqn)
         end do
#endif
      end subroutine HexMesh_UpdateMPIFacesAviscFlux
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_GatherMPIFacesSolution(self, nEqn)
         implicit none
         class(HexMesh)    :: self
         integer, intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter
         integer, parameter :: otherSide(2) = (/2,1/)

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************
!        Gather solution
!        ***************
!
         do domain = 1, MPI_Process % nProcs
!
!           **************************************
!           Wait until messages have been received
!           **************************************
!
            call self % MPIfaces % faces(domain) % WaitForSolution

            counter = 1
            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  f % storage(otherSide(thisSide)) % Q(:,i,j) = self % MPIfaces % faces(domain) % Qrecv(counter:counter+nEqn-1)
                  counter = counter + nEqn
               end do               ; end do
               end associate
            end do
         end do
#endif
      end subroutine HexMesh_GatherMPIFacesSolution

      subroutine HexMesh_GatherMPIFacesGradients(self, nEqn)
         implicit none
         class(HexMesh)      :: self
         integer, intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter
         integer, parameter :: otherSide(2) = (/2,1/)

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************
!        Gather solution
!        ***************
!
         do domain = 1, MPI_Process % nProcs

!
!           **************************************
!           Wait until messages have been received
!           **************************************
!
            call self % MPIfaces % faces(domain) % WaitForGradients

            counter = 1
            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  f % storage(otherSide(thisSide)) % U_x(:,i,j) = self % MPIfaces % faces(domain) % U_xyzrecv(counter:counter+nEqn-1)
                  counter = counter + nEqn
               end do               ; end do

               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  f % storage(otherSide(thisSide)) % U_y(:,i,j) = self % MPIfaces % faces(domain) % U_xyzrecv(counter:counter+nEqn-1)
                  counter = counter + nEqn
               end do               ; end do

               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  f % storage(otherSide(thisSide)) % U_z(:,i,j) = self % MPIfaces % faces(domain) % U_xyzrecv(counter:counter+nEqn-1)
                  counter = counter + nEqn
               end do               ; end do
               end associate
            end do
         end do
#endif
      end subroutine HexMesh_GatherMPIFacesGradients
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_GatherMPIFacesAviscflux(self, nEqn)
         implicit none
         class(HexMesh)      :: self
         integer, intent(in) :: nEqn
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: mpifID, fID, thisSide, domain
         integer            :: i, j, counter
         integer, parameter :: otherSide(2) = [2, 1]

         if ( .not. MPI_Process % doMPIAction ) return
!
!        ***************
!        Gather solution
!        ***************
!
         do domain = 1, MPI_Process % nProcs
!
!           **************************************
!           Wait until messages have been received
!           **************************************
!
            call self % MPIfaces % faces(domain) % WaitForAviscflux

            counter = 1
            do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
               fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
               thisSide = self % MPIfaces % faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  f % storage(otherSide(thisSide)) % Aviscflux(:,i,j) = &
                     self % MPIfaces % faces(domain) % AviscfluxRecv(counter:counter+nEqn-1)
                  counter = counter + nEqn
               end do               ; end do
               end associate
            end do
         end do
#endif
      end subroutine HexMesh_GatherMPIFacesAviscflux
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_PrepareForIO(self)
         implicit none
         class(HexMesh)    :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer              :: ierr, eID
         integer, allocatable :: elementSizes(:)
         integer, allocatable :: allElementSizes(:)
         integer, allocatable :: allElementsOffset(:)
!
!        Get each element size
!        ---------------------
         allocate(elementSizes(self % no_of_allElements))

         elementSizes = 0     ! Default value 0 to use allreduce with SUM
         do eID = 1, self % no_of_elements
            associate(e => self % elements(eID))
            elementSizes(e % globID) = product( e % Nxyz + 1)
            end associate
         end do
!
!        Gather the rest of the mesh values if it is a partition
!        -------------------------------------------------------
         allocate(allElementSizes(self % no_of_allElements))

         allElementSizes = 0
         if ( (MPI_Process % doMPIAction) ) then
#ifdef _HAS_MPI_
            call mpi_allreduce(elementSizes, allElementSizes, &
                               self % no_of_allElements, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
         else
            allElementSizes = elementSizes

         end if
!
!        Get all elements offset: the accumulation of allElementSizes
!        -----------------------
         allocate(allElementsOffset(self % no_of_allElements))

         allElementsOffset(1) = 0
         do eID = 2, self % no_of_allElements
            allElementsOffset(eID) = allElementsOffset(eID-1) + allElementSizes(eID-1)
         end do
!
!        Assign the results to the elements
!        ----------------------------------
         do eID = 1, self % no_of_elements
            associate(e => self % elements(eID))
            e % offsetIO = allElementsOffset(e % globID)
            end associate
         end do
!
!        Free memory
!        -----------
         deallocate(elementSizes, allElementSizes, allElementsOffset)

      end subroutine HexMesh_PrepareForIO
!
!////////////////////////////////////////////////////////////////////////
!
!     --------------------------------------------------------------------
!     This subroutine describes the loaded mesh
!     --------------------------------------------------------------------
      SUBROUTINE HexMesh_Describe( self , fileName, bFaceOrder )
      USE Headers
      IMPLICIT NONE
      !-arguments------------------------------------------
      CLASS(HexMesh)      :: self
      CHARACTER(LEN=*)    :: fileName
      integer, intent(in) :: bFaceOrder
      !-local-variables------------------------------------
      integer           :: ierr
      integer           :: zoneID
      integer           :: no_of_bdry_faces
      integer           :: no_of_faces
      integer, allocatable :: facesPerZone(:)
      character(len=LINE_LENGTH) :: str
      !----------------------------------------------------

      allocate ( facesPerZone(size(self % zones)) )

!     Gather information
!     ------------------

      if (  MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
         do zoneID = 1, size(self % zones)
            call mpi_reduce ( self % zones(zoneID) % no_of_faces, facesPerZone(zoneID) , 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
         end do

         no_of_bdry_faces = sum(facesPerZone)
         no_of_faces      = (6*self % no_of_allElements + no_of_bdry_faces)/2
#endif
      else
         do zoneID = 1, size(self % zones)
            facesPerZone(zoneID) = self % zones(zoneID) % no_of_faces
         end do

         no_of_bdry_faces = sum(facesPerZone)
         no_of_faces = size ( self % faces )
      end if


!     Describe the mesh
!     -----------------

      if ( .not. MPI_Process % isRoot ) return

      write(STD_OUT,'(/)')
      call Section_Header("Mesh information")
      write(STD_OUT,'(/)')

      call SubSection_Header('Mesh file "' // trim(fileName) // '".')

      write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of elements: " , self % no_of_allElements
      write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of faces: " , no_of_faces

      write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of boundary faces: " , no_of_bdry_faces
      write(STD_OUT,'(30X,A,A28,I10)') "->" , "Order of curved faces: " , bFaceOrder
      write(STD_OUT,'(30X,A,A28,L10)') "->" , "2D extruded mesh: " , self % meshIs2D

!     Describe the zones
!     ------------------
      write(STD_OUT,'(/)')
      call Section_Header("Creating zones")
      write(STD_OUT,'(/)')

      do zoneID = 1, size(self % zones)
         write(str,'(A,I0,A,A)') "Zone ", zoneID, " for boundary: ",trim(self % zones(zoneID) % Name)
         call SubSection_Header(trim(str))
         write(STD_OUT,'(30X,A,A28,I0)') "->", ' Number of faces: ', facesPerZone(zoneID)
         call BCs(zoneID) % bc % Describe
         write(STD_OUT,'(/)')
      end do

!     Finish up
!     ---------
      deallocate ( facesPerZone )

      END SUBROUTINE HexMesh_Describe
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DescribeMeshPartition( self )
      USE Headers
      use PartitionedMeshClass
      IMPLICIT NONE
!
!--------------------------------------------------------------------
!  This subroutine describes the loaded mesh partition
!--------------------------------------------------------------------
!
!
!     ------------------
!     External variables
!     ------------------
!
      CLASS(HexMesh)    :: self
#ifdef _HAS_MPI_
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER           :: fID, zoneID, rank, ierr
      INTEGER           :: no_of_bdryfaces, no_of_mpifaces
      integer           :: no_of_elementsP(MPI_Process % nProcs)
      integer           :: no_of_facesP(MPI_Process % nProcs)
      integer           :: no_of_bfacesP(MPI_Process % nProcs)
      integer           :: no_of_mpifacesP(MPI_Process % nProcs)
      character(len=64) :: partitionID

      if ( .not. MPI_Process % doMPIAction ) return

      no_of_bdryfaces = 0
      no_of_mpifaces  = 0

      do fID = 1 , size ( self % faces )
         if ( self % faces(fID) % faceType .eq. HMESH_BOUNDARY) then
            no_of_bdryfaces = no_of_bdryfaces + 1
         elseif ( self % faces(fID) % faceType .eq. HMESH_MPI ) then
            no_of_mpifaces = no_of_mpifaces + 1
         end if
      end do
!
!     Share all quantities to the root process
!     ----------------------------------------
      call mpi_gather(size(self % elements) , 1 , MPI_INT , no_of_elementsP , 1 , MPI_INT , 0 , MPI_COMM_WORLD , ierr)
      call mpi_gather(size(self % faces)    , 1 , MPI_INT , no_of_facesP    , 1 , MPI_INT , 0 , MPI_COMM_WORLD , ierr)
      call mpi_gather(no_of_bdryfaces       , 1 , MPI_INT , no_of_bfacesP   , 1 , MPI_INT , 0 , MPI_COMM_WORLD , ierr)
      call mpi_gather(no_of_mpifaces        , 1 , MPI_INT , no_of_mpifacesP , 1 , MPI_INT , 0 , MPI_COMM_WORLD , ierr)

      if ( .not. MPI_Process % isRoot ) return

      write(STD_OUT,'(/)')
      call Section_Header("Mesh partitions")
      write(STD_OUT,'(/)')

      write(STD_OUT,'(10X,A21)', advance='no') "Partitioning method: "

      select case (MPI_Partitioning)
         case (METIS_PARTITIONING)
            write(STD_OUT,'(A)') 'METIS'
         case (SFC_PARTITIONING)
            write(STD_OUT,'(A)') 'Space-filling curve'
      end select
      write(STD_OUT,*)

      do rank = 1, MPI_Process % nProcs

         write(partitionID,'(A,I0)') "Partition ", rank
         call SubSection_Header(trim(partitionID))

         write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of elements: " , no_of_elementsP(rank)
         write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of faces: " , no_of_facesP(rank)
         write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of boundary faces: " , no_of_bfacesP(rank)
         write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of mpi faces: " , no_of_mpifacesP(rank)

      end do
#endif

      END SUBROUTINE DescribeMeshPartition

!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE WriteCoordFile(self,nEqn, FileName)
         USE PhysicsStorage
         IMPLICIT NONE
!
!        -----------------------------------------------------------------
!        This subroutine writes a *.coo file containing all the mesh nodes
!        that can be used for eigenvalue analysis using the TAUev code
!        -----------------------------------------------------------------
!
         !--------------------------------------------------------
         CLASS(HexMesh)       :: self        !<  this mesh
         integer              :: nEqn
         CHARACTER(len=*)     :: FileName    !<  ...
         !--------------------------------------------------------
         INTEGER              :: NumOfElem
         INTEGER              :: i, j, k, el, Nx, Ny, Nz, ndof, cooh
         !--------------------------------------------------------

         NumOfElem = SIZE(self % elements)
!
!        ------------------------------------------------------------------------
!        Determine the number of degrees of freedom
!        ------------------------------------------------------------------------
!
         ndof = 0
         DO el = 1, NumOfElem
            Nx = self % elements(el) % Nxyz(1)
            Ny = self % elements(el) % Nxyz(2)
            Nz = self % elements(el) % Nxyz(3)
            ndof = ndof + (Nx + 1)*(Ny + 1)*(Nz + 1)*nEqn
         END DO

         OPEN(newunit=cooh, file=FileName, action='WRITE')

         WRITE(cooh,*) ndof, ndim   ! defined in PhysicsStorage
         DO el = 1, NumOfElem
            Nx = self % elements(el) % Nxyz(1)
            Ny = self % elements(el) % Nxyz(2)
            Nz = self % elements(el) % Nxyz(3)
            DO k = 0, Nz
               DO j = 0, Ny
                  DO i = 0, Nx
                     WRITE(cooh,*) self % elements(el) % geom % x(1,i,j,k), &
                                   self % elements(el) % geom % x(2,i,j,k), &
                                   self % elements(el) % geom % x(3,i,j,k)
                  END DO
               END DO
            END DO
         END DO

         CLOSE(cooh)

      END SUBROUTINE WriteCoordFile
!
!//////////////////////////////////////////////////////////////////////////////
!
!     This subroutine checks if the mesh is a 2D extruded mesh and gets the 2D direction
!     in the local frame for each element
!     --------------------------------------------------------------------
!
!//////////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_CheckIfMeshIs2D(self)
         implicit none
         !-arguments---------------------------------------
         class(HexMesh),   intent(inout) :: self
         !-local-variables---------------------------------
         integer  :: eID, nID, no_of_orientedNodes
         integer  :: dir
         integer  :: ierr                             ! Error for MPI calls
         integer  :: face1Nodes(NODES_PER_FACE)
         integer  :: face2Nodes(NODES_PER_FACE)
         integer  :: no_of_orientedElems(NDIM)
         logical  :: meshExtrudedIn(NDIM)
         logical  :: meshExtrudedInLocal(NDIM)
         real(kind=RP)  :: xNodesF1(NDIM,NODES_PER_FACE)
         real(kind=RP)  :: xNodesF2(NDIM,NODES_PER_FACE)
         real(kind=RP)  :: dx(NDIM,NODES_PER_FACE)
         real(kind=RP), parameter   :: d2D(NDIM,NDIM) = reshape(  (/1._RP, 0._RP, 0._RP, &
                                                                    0._RP, 1._RP, 0._RP, &
                                                                    0._RP, 0._RP, 1._RP /) , (/3,3/) )
         !-------------------------------------------------

         no_of_orientedElems = 0

         elem_loop: do eID = 1, self % no_of_elements
            associate(e => self % elements(eID))
!
!           *****************************************
!           Check if the direction is xi (Left,Right)
!           *****************************************
!
!           Get both face node IDs
!           ----------------------
            face1Nodes = e % nodeIDs(localFaceNode(:, ELEFT))
            face2Nodes = e % nodeIDs(localFaceNode(:, ERIGHT))
!
!           Get the nodes coordinates
!           -------------------------
            do nID = 1, NODES_PER_FACE
               xNodesF1(:,nID) = self % nodes(face1Nodes(nID)) % x
               xNodesF2(:,nID) = self % nodes(face2Nodes(nID)) % x

            end do
!
!           Compute the delta x vectors
!           ---------------------------
            dx = xNodesF2 - xNodesF1
!
!           Normalize!
!           ----------
            do nID = 1, NODES_PER_FACE
               dx(:,nID) = dx(:,nID) / norm2( dx(:,nID) )
            end do
!
!           Check how many delta x vectors are parallel to the the x, y or z axis
!           ---------------------------------------------------------------

            do dir = 1, NDIM
               no_of_orientedNodes = 0
               do nID = 1, NODES_PER_FACE
                  if ( almostEqual(abs(dot_product(dx(:,nID),d2D(:,dir))),1.0_RP) ) then
                     no_of_orientedNodes = no_of_orientedNodes + 1
                  end if
               end do

               if ( no_of_orientedNodes .eq. 4 ) then
!
!                 This is (at least one of) the 2D direction(s)
!                 ---------------------------------------------
                  e % dir2D = IX
                  no_of_orientedElems(dir) = no_of_orientedElems(dir) + 1
                  e % globDir(dir) = IX

               end if
            end do

!
!           *****************************************
!           Check if the direction is eta (Front,Back)
!           *****************************************
!
!           Get both face node IDs
!           ----------------------
            face1Nodes = e % nodeIDs(localFaceNode(:, EFRONT))
            face2Nodes = e % nodeIDs(localFaceNode(:, EBACK))
!
!           Get the nodes coordinates
!           -------------------------
            do nID = 1, NODES_PER_FACE
               xNodesF1(:,nID) = self % nodes(face1Nodes(nID)) % x
               xNodesF2(:,nID) = self % nodes(face2Nodes(nID)) % x
            end do
!
!           Compute the delta x vectors
!           ---------------------------
            dx = xNodesF2 - xNodesF1
!
!           Normalize!
!           ----------
            do nID = 1, NODES_PER_FACE
               dx(:,nID) = dx(:,nID) / norm2( dx(:,nID) )
            end do
!
!           Check how many delta x vectors are parallel to the 2D direction
!           ---------------------------------------------------------------
            do dir = 1, NDIM
               no_of_orientedNodes = 0
               do nID = 1, NODES_PER_FACE
                  if ( almostEqual(abs(dot_product(dx(:,nID),d2D(:,dir))),1.0_RP) ) then
                     no_of_orientedNodes = no_of_orientedNodes + 1
                  end if
               end do

               if ( no_of_orientedNodes .eq. 4 ) then
!
!                 This is (at least one of) the 2D direction(s)
!                 ---------------------------------------------
                  if (e % dir2D == IX) then
                     e % dir2D = IXY
                  else
                     e % dir2D = IY
                  end if
                  no_of_orientedElems(dir) = no_of_orientedElems(dir) + 1
                  e % globDir(dir) = IY

               end if

            end do
!
!           *****************************************
!           Check if the direction is zeta (Bottom,Top)
!           *****************************************
!
!           Get both face node IDs
!           ----------------------
            face1Nodes = e % nodeIDs(localFaceNode(:, EBOTTOM))
            face2Nodes = e % nodeIDs(localFaceNode(:, ETOP))
!
!           Get the nodes coordinates
!           -------------------------
            do nID = 1, NODES_PER_FACE
               xNodesF1(:,nID) = self % nodes(face1Nodes(nID)) % x
               xNodesF2(:,nID) = self % nodes(face2Nodes(nID)) % x
            end do
!
!           Compute the delta x vectors
!           ---------------------------
            dx = xNodesF2 - xNodesF1
!
!           Normalize!
!           ----------
            do nID = 1, NODES_PER_FACE
               dx(:,nID) = dx(:,nID) / norm2( dx(:,nID) )
            end do
!
!           Check how many delta x vectors are parallel to the 2D direction
!           ---------------------------------------------------------------
            do dir = 1, NDIM
               no_of_orientedNodes = 0
               do nID = 1, NODES_PER_FACE
                  if ( almostEqual(abs(dot_product(dx(:,nID),d2D(:,dir))),1.0_RP) ) then
                     no_of_orientedNodes = no_of_orientedNodes + 1
                  end if
               end do

               if ( no_of_orientedNodes .eq. 4 ) then
!
!                 This is (at least one of) the 2D direction(s)
!                 ---------------------------------------------
                  select case (e % dir2D)
                     case (IX)
                        e % dir2D = IXZ
                     case (IY)
                        e % dir2D = IYZ
                     case (IXY)
                        e % dir2D = IXYZ
                     case default
                        e % dir2D = IZ
                  end select

                  no_of_orientedElems(dir) = no_of_orientedElems(dir) + 1
                  e % globDir(dir) = IZ

               end if
            end do

            end associate
         end do elem_loop

         meshExtrudedIn = ( no_of_orientedElems == self % no_of_elements )

!        MPI communication
!        -----------------
#if _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
            meshExtrudedInLocal = meshExtrudedIn
            call mpi_allreduce ( meshExtrudedInLocal, meshExtrudedIn, NDIM, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr )
         end if
#endif

         if ( any(meshExtrudedIn) ) then
            self % meshIs2D = .TRUE.

            if ( all(meshExtrudedIn) ) then
               self % dir2D = IXYZ
            elseif ( meshExtrudedIn(IX) .and. meshExtrudedIn(IY) ) then
               self % dir2D = IXY
            elseif ( meshExtrudedIn(IX) .and. meshExtrudedIn(IZ) ) then
               self % dir2D = IXZ
            elseif ( meshExtrudedIn(IY) .and. meshExtrudedIn(IZ) ) then
               self % dir2D = IYZ
            elseif ( meshExtrudedIn(IX) ) then
               self % dir2D = IX
            elseif ( meshExtrudedIn(IY) ) then
               self % dir2D = IY
            elseif ( meshExtrudedIn(IZ) ) then
               self % dir2D = IZ
            end if
         end if

      end subroutine HexMesh_CheckIfMeshIs2D
!
!//////////////////////////////////////////////////////////////////////////////
!
!     If the mesh is a 2D extruded mesh, this subroutine sets the polynomial order
!     to zero in the corresponding direction (HexMesh_CheckIfMeshIs2D must have been called beforehand)
!     --------------------------------------------------------------------
!
!//////////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_CorrectOrderFor2DMesh(self, dir2D,order_)
         implicit none
         class(HexMesh),   intent(inout) :: self
         integer,          intent(in)    :: dir2D
         integer, intent(in), optional   :: order_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: eID, nID, no_of_orientedNodes, order
         integer  :: face1Nodes(NODES_PER_FACE)
         integer  :: face2Nodes(NODES_PER_FACE)
         logical  :: rightDir
         real(kind=RP)  :: d2D(NDIM)
         real(kind=RP)  :: xNodesF1(NDIM,NODES_PER_FACE)
         real(kind=RP)  :: xNodesF2(NDIM,NODES_PER_FACE)
         real(kind=RP)  :: dx(NDIM,NODES_PER_FACE)

         if ( present(order_) ) then
            order = order_
         else
            order = 0
         end if

         if (self % meshIs2D) then
            select case (dir2D)
               case (IX)
                  rightDir = any (self % dir2D == [IX, IXY, IXZ, IXYZ] )
               case (IY)
                  rightDir = any (self % dir2D == [IY, IXY, IYZ, IXYZ] )
               case (IZ)
                  rightDir = any (self % dir2D == [IZ, IXZ, IYZ, IXYZ] )
            end select
            if (.not. rightDir) then
               print*, "The mesh does not seem to be 2D for the selected direction"
               errorMessage(STD_OUT)
               return
            end if
         else
            print*, "The mesh does not seem to be 2D"
            errorMessage(STD_OUT)
         end if

         do eID = 1, self % no_of_elements
            associate(e => self % elements(eID))

            select case (e % globDir(dir2D))
               case (IX)
                  e % Nxyz(1) = order
                  self % Nx(eID) = order
               case (IY)
                  e % Nxyz(2) = order
                  self % Ny(eID) = order
               case (IZ)
                  e % Nxyz(3) = order
                  self % Nz(eID) = order
            end select

            end associate
         end do

      end subroutine HexMesh_CorrectOrderFor2DMesh
!
!////////////////////////////////////////////////////////////////////////
!
!        Set element connectivities
!        --------------------------
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_SetConnectivitiesAndLinkFaces(self,nodes,facesList)
         implicit none
         !-arguments----------------------------------------------------
         class(HexMesh), target, intent(inout) :: self
         integer               , intent(in)    :: nodes
         integer, optional     , intent(in)    :: facesList(:)
         !-local-variables----------------------------------------------
         integer  :: fID, SideL, SideR, side, counter
         integer  :: NelL(2), NelR(2)    ! Polynomial orders on left and right of a face
         integer  :: Nel(2,2)            ! Polynomial orders on left and right of a face - for MPI
         integer  :: globID
         integer  :: Nxyz(NDIM)
         integer  :: domain, MPI_NDOFS(MPI_Process % nProcs), mpifID
         integer  :: num_of_Faces, ii
         integer, parameter :: other(2) = [2, 1]
         !--------------------------------------------------------------

         if ( present(facesList) ) then
            num_of_Faces = size(facesList)
         else
            num_of_Faces = size(self % faces)
         end if
!
!        ---------------------------------------------------------------------------
!        Send polynomial orders across MPI faces
!           -> Currently doing for all MPI faces (not taking into account facesList)
!        ---------------------------------------------------------------------------
!
         if (mpi_partition % Constructed) call self % UpdateMPIFacesPolynomial
!
!        ------------------------
!        Link faces with elements
!        ------------------------
!
         do ii = 1, num_of_Faces

            if ( present(facesList) ) then
               fID = facesList(ii)
            else
               fID = ii
            end if
            associate (  f => self % faces(fID)   )

            select case (f % faceType)
            case (HMESH_INTERIOR)
               associate(eL => self % elements(f % elementIDs(1)), &
                         eR => self % elements(f % elementIDs(2))   )
!
!              Get polynomial orders of elements
!              ---------------------------------
               NelL = eL % Nxyz(axisMap(:, f % elementSide(1)))
               NelR = eR % Nxyz(axisMap(:, f % elementSide(2)))
!
!              Fill connectivity of element type
!              ---------------------------------
               SideL = f % elementSide(1)
               SideR = f % elementSide(2)
!
!              Construct connectivity
!              ----------------------
               eL % NumberOfConnections(SideL) = 1
               call eL % Connection(SideL) % Construct(eR % GlobID, eR % Nxyz)

               eR % NumberOfConnections(SideR) = 1
               call eR % Connection(SideR) % Construct(eL % GlobID, eL % Nxyz)
               end associate

               call f % LinkWithElements(NelL, NelR, nodes)

            case (HMESH_BOUNDARY)
               associate(eL => self % elements(f % elementIDs(1)))
!
!              Get polynomial orders of elements
!              ---------------------------------
               NelL = eL % Nxyz(axisMap(:, f % elementSide(1)))
               NelR = NelL

               ! Default NumberOfConnections = 0
               end associate

               call f % LinkWithElements(NelL, NelR, nodes)

!           case (HMESH_MPI): Do nothing. MPI faces are constructed in the next step.

            end select


            end associate
         end do
!
!        -----------------------------------------------
!        Gather faces polynomial and link MPI faces
!        -> All, not only the faces included in faceList
!        -----------------------------------------------
!
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction .and. mpi_partition % Constructed  )  then

            do domain = 1, MPI_Process % nProcs
!
!              Wait until messages have been received
!              --------------------------------------
!
               call self % MPIfaces % faces(domain) % WaitForN

               counter = 1
               do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
                  fID  = self % MPIfaces % faces(domain) % faceIDs(mpifID)
                  side = self % MPIfaces % faces(domain) % elementSide(mpifID)   ! face side 1/2

                  associate( f => self % faces(fID) )
                  associate( e => self % elements(maxval(f % elementIDs)) )

                  sideL = f % elementSide(side)                                  ! element side 1/2/3/4/5/6

                  Nel(:,      side ) = e % Nxyz(axisMap(:,sideL))
                  Nel(:,other(side)) = self % MPIfaces % faces(domain) % Nrecv(counter:counter+1)

                  call f % LinkWithElements(Nel(:,1), Nel(:,2), nodes)

                  Nxyz   = self % MPIfaces % faces(domain) % Nrecv(counter+2:counter+4)
                  globID = self % MPIfaces % faces(domain) % Nrecv(counter+5)

                  e % NumberOfConnections (sideL) = 1
                  call e % Connection(sideL) % construct (globID,Nxyz)
                  counter = counter + 6

                  end associate
                  end associate
               end do
            end do
         end if
#endif
!
!        --------------------------
!        Allocate MPI Faces storage
!        --------------------------
!
         if ( MPI_Process % doMPIAction ) then
            if ( .not. allocated(self % MPIfaces % faces) ) return
#if _HAS_MPI_
            MPI_NDOFS = 0

            do domain = 1, MPI_Process % nProcs
               do mpifID = 1, self % MPIfaces % faces(domain) % no_of_faces
                  fID = self % MPIfaces % faces(domain) % faceIDs(mpifID)
                  associate( fc => self % faces(fID) )
                  MPI_NDOFS(domain) = MPI_NDOFS(domain) + product(fc % Nf + 1)
                  end associate
               end do
            end do

#if defined(NAVIERSTOKES)
            call ConstructMPIFacesStorage(self % MPIfaces, NCONS, NGRAD, MPI_NDOFS)
#elif defined(INCNS)
            call ConstructMPIFacesStorage(self % MPIfaces, NCONS, NCONS, MPI_NDOFS)
#elif defined(CAHNHILLIARD)
            call ConstructMPIFacesStorage(self % MPIfaces, NCOMP, NCOMP, MPI_NDOFS)
#endif

#endif
         end if

      end subroutine HexMesh_SetConnectivitiesAndLinkFaces
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_UpdateFacesWithPartition(self, partition, nAllElems, global2LocalIDs)
!
!        **************************************************
!        This subroutine casts HMESH_UNDEFINED faces which
!        are a partition boundary face as HMESH_MPI
!        **************************************************
!
         use PartitionedMeshClass
         implicit none
         class(HexMesh)    :: self
         class(PartitionedMesh_t),  intent(in)  :: partition
         integer,                   intent(in)  :: nAllElems
         integer,                   intent(in)  :: global2LocalIDs(nAllElems)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: k, eID, bFace, side, eSide, fID, domain
         integer  :: no_of_mpifaces(MPI_Process % nProcs)
         integer, parameter  :: otherSide(2) = (/2,1/)
         integer, parameter  :: invRot(1:4,0:7) = reshape( (/ 1, 2, 3, 4, &
                                                              4, 1, 2, 3, &
                                                              3, 4, 1, 2, &
                                                              2, 3, 4, 1, &
                                                              1, 4, 3, 2, &
                                                              2, 1, 4, 3, &
                                                              3, 2, 1, 4, &
                                                              4, 3, 2, 1 /), (/4,8/) )
!
!        First get how many faces are shared with each other partition
!        -------------------------------------------------------------
         no_of_mpifaces = 0
         do bFace = 1, partition % no_of_mpifaces
            domain = partition % mpiface_sharedDomain(bFace)
            no_of_mpifaces(domain) = no_of_mpifaces(domain) + 1
         end do
!
!        ---------------
!        Allocate memory
!        ---------------
!
         do domain = 1, MPI_Process % nProcs
            if ( no_of_mpifaces(domain) .ne. 0 ) then
               call self % MPIfaces % faces(domain) % Construct(no_of_mpifaces(domain))
            end if
         end do
!
!        -------------
!        Assign values
!        -------------
!
         no_of_mpifaces = 0
         do bFace = 1, partition % no_of_mpifaces
!
!           Gather the face, and the relative position w.r.t. its element
!           -------------------------------------------------------------
            eID = global2LocalIDs(partition % mpiface_elements(bFace))
            side = partition % element_mpifaceSide(bFace)
            eSide = partition % mpiface_elementSide(bFace)
            fID = self % elements(eID) % faceIDs(side)
!
!           Change the face to a HMESH_MPI
!           ------------------------------
            associate(f => self % faces(fID))
            f % faceType = HMESH_MPI
            f % rotation = partition % mpiface_rotation(bFace)
            f % elementIDs(eSide) = eID
            f % elementIDs(otherSide(eSide)) = HMESH_NONE   ! This makes sense since elementIDs are in local numbering...
            f % elementSide(eSide) = side
            f % elementSide(otherSide(eSide)) = partition % element_mpifaceSideOther(bFace)

            if (eSide == RIGHT) f % nodeIDs = f % nodeIDs (invRot (:,f % rotation) )
            end associate

            self % elements(eID) % faceSide(side) = eSide

!
!           Create MPI Face
!           ---------------
            domain = partition % mpiface_sharedDomain(bFace)
            no_of_mpifaces(domain) = no_of_mpifaces(domain) + 1
            self % MPIfaces % faces(domain) % faceIDs(no_of_mpifaces(domain)) = fID
            self % MPIfaces % faces(domain) % elementSide(no_of_mpifaces(domain)) = eSide

         end do

      end subroutine HexMesh_UpdateFacesWithPartition
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------------------
!     Construct geometry of faces and elements
!     -> This routine guarantees that the mapping is subparametric or isoparametric
!     -----------------------------------------------------------------------------
      subroutine HexMesh_ConstructGeometry(self, facesList, elementList)
         implicit none
         !--------------------------------
         class(HexMesh)    , intent(inout) :: self
         integer, optional , intent(in)    :: facesList(:)
         integer, optional , intent(in)    :: elementList(:)
         !--------------------------------
         integer                       :: num_of_elements, num_of_faces
         integer                       :: i, ii
         integer                       :: fID, eID                        ! Face and element counters
         integer                       :: eIDLeft, eIDRight, e            ! Element IDs on left and right of a face
         integer                       :: SideIDL, SideIDR, side          ! Side of elements on left and right of a face
         integer                       :: buffer                          ! A temporal variable
         integer                       :: NSurfL(2), NSurfR(2), NSurf(2)  ! Polynomial order the face was constructef with
         integer                       :: Nelf(2), Nel(2), rot
         integer                       :: CLN(2)                          ! Chebyshev-Lobatto face orders
         REAL(KIND=RP)                 :: corners(NDIM,NODES_PER_ELEMENT) ! Variable just for initializing purposes
         real(kind=RP)   , allocatable :: faceCL(:,:,:)                   ! Coordinates of the Chebyshev-Lobatto nodes on the face
         type(SurfInfo_t), allocatable :: SurfInfo(:)                     ! Local copy of surf info that can be modified
         type(TransfiniteHexMap), pointer :: hexMap, hex8Map, genHexMap
         !--------------------------------
         logical                    :: isConforming   ! Is the representation conforming on a boundary?
         integer                    :: zoneID, zonefID, nZones
         integer, allocatable       :: bfOrder_local(:)  ! Polynomial order on a zone (partition wise)
         integer, allocatable       :: bfOrder(:)        ! Polynomial order on a zone (global)
         integer                    :: ierr              ! Error for MPI calls
         integer                    :: thisSide          ! The side of the MPI face that corresponds to current partition
         !--------------------------------

         corners = 0._RP

         if ( present(elementList) ) then
            num_of_elements = size(elementList)
         else
            num_of_elements = size(self % elements)
         end if
         if ( present(facesList) ) then
            num_of_faces = size(facesList)
         else
            num_of_faces = size(self % faces)
         end if

!
!        Generate a local copy of SurfInfo
!        ---------------------------------
         allocate ( SurfInfo(self % no_of_elements) )
         do ii=1, num_of_elements
            if ( present(elementList) ) then
               eID = elementList(ii)
            else
               eID = ii
            end if

            SurfInfo(eID) = self % elements(eID) % SurfInfo
         end do
!
!        ***********************************************************************
!        Find the polynomial order of the boundaries for 3D nonconforming meshes
!        -> In 3D meshes we force all the faces on a boundary to be mapped with the
!           same polynomial order, in order to have "water-tight" meshes. This condition
!           can be sometimes too strict...
!        -> When the representation is p-nonconforming on a boundary, the mapping
!           must be of order P <= min(N)/2. This is the general sufficient condition.
!        *************************************************************1**********
!
         if (self % anisotropic .and. (.not. self % meshIs2D) ) then
            nZones = size(self % zones)
            allocate ( bfOrder (nZones) )
            bfOrder = huge(bfOrder) ! Initialize to a big number

            do zoneID=1, nZones
               if (self % zones(zoneID) % no_of_faces == 0 ) cycle

!              Get the minimum polynomial order in this zone
!              ---------------------------------------------
               do zonefID = 1, self % zones(zoneID) % no_of_faces
                  fID = self % zones(zoneID) % faces(zonefID)

                  associate( f => self % faces(fID) )
                  bfOrder(zoneID) = min(bfOrder(zoneID),f % NfLeft(1),f % NfLeft(2))
                  end associate
               end do

!           MPI communication
!           -----------------
#if _HAS_MPI_
            end do

            if ( MPI_Process % doMPIAction ) then
               allocate ( bfOrder_local(nZones) )
               bfOrder_local = bfOrder
               call mpi_allreduce ( bfOrder_local, bfOrder, nZones, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr )
               deallocate ( bfOrder_local )
            end if

            do zoneID=1, nZones
#endif

!              Select the BC poynomial order
!              -----------------------------

               if ( self % ConformingOnZone(zoneID) .or. self % ignoreBCnonConformities) then
                  !bfOrder(zoneID) = bfOrder(zoneID)
               else
                  bfOrder(zoneID) = bfOrder(zoneID)/2
                  if ( bfOrder(zoneID) < 1 ) then
                     write(STD_OUT,*) 'ERROR :: The chosen polynomial orders are too low to represent the boundaries accurately'
                     write(STD_OUT,*) '      :: Nonconforming representations on boundaries need N>=2'
                     stop
                  end if
               end if

               call NodalStorage (bfOrder(zoneID)) % construct (self % nodeType, bfOrder(zoneID))

            end do
         end if
!
!        **************************************************************
!        Check surfaces' integrity and adapt them to the solution order
!        **************************************************************
!
         do ii=1, num_of_faces
            if ( present(facesList) ) then
               fID = facesList(ii)
            else
               fID = ii
            end if
            associate( f => self % faces(fID) )

!
!           Check if the surfaces description in mesh file is consistent
!           ------------------------------------------------------------
            select case (f % faceType)

            case (HMESH_INTERIOR)
               eIDLeft  = f % elementIDs(1)
               SideIDL  = f % elementSide(1)
               NSurfL   = SurfInfo(eIDLeft) % facePatches(SideIDL) % noOfKnots - 1

               eIDRight = f % elementIDs(2)
               SideIDR  = f % elementSide(2)
               NSurfR   = SurfInfo(eIDRight) % facePatches(SideIDR) % noOfKnots - 1

!              If both surfaces are of order 1.. There's no need to continue analyzing face
!              ----------------------------------------------------------------------------
               if     ((SurfInfo(eIDLeft)  % IsHex8) .and. (SurfInfo(eIDRight) % IsHex8)) then
                  cycle
               elseif ((SurfInfo(eIDLeft)  % IsHex8) .and. all(NSurfR == 1) ) then
                  cycle
               elseif ((SurfInfo(eIDRight) % IsHex8) .and. all(NSurfL == 1) ) then
                  cycle
               elseif (all(NSurfL == 1) .and. all(NSurfR == 1) ) then
                  cycle
               elseif (any(NSurfL /= NSurfR)) then ! Only works for mesh files with isotropic boundary orders
                  write(STD_OUT,*) 'WARNING: Curved face definitions in mesh are not consistent.'
                  write(STD_OUT,*) '   Face:    ', fID
                  write(STD_OUT,*) '   Elements:', f % elementIDs
                  write(STD_OUT,*) '   N Left:  ', SurfInfo(eIDLeft)  % facePatches(SideIDL) % noOfKnots - 1
                  write(STD_OUT,*) '   N Right: ', SurfInfo(eIDRight) % facePatches(SideIDR) % noOfKnots - 1
               end if

               CLN(1) = min(f % NfLeft(1),f % NfRight(1))
               CLN(2) = min(f % NfLeft(2),f % NfRight(2))
!
!              Adapt the curved face order to the polynomial order
!              ---------------------------------------------------
               if ( any(CLN < NSurfL) ) then
                  allocate(faceCL(1:3,CLN(1)+1,CLN(2)+1))
                  call ProjectFaceToNewPoints(SurfInfo(eIDLeft) % facePatches(SideIDL), CLN(1), NodalStorage(CLN(1)) % xCGL, &
                                                                                        CLN(2), NodalStorage(CLN(2)) % xCGL, faceCL)
                  call SurfInfo(eIDLeft) % facePatches(SideIDL) % Destruct()
                  call SurfInfo(eIDLeft) % facePatches(SideIDL) % Construct(NodalStorage(CLN(1)) % xCGL, &
                                                                            NodalStorage(CLN(2)) % xCGL,faceCL)
                  deallocate(faceCL)
               end if

               select case ( f % rotation )
               case ( 1, 3, 4, 6 ) ! Local x and y axis are perpendicular
                  if (CLN(1) /= CLN(2)) then
                     buffer = CLN(1)
                     CLN(1) = CLN(2)
                     CLN(2) = buffer
                  end if
               end select

               if ( any(CLN < NSurfR) ) then      
                  allocate(faceCL(1:3,CLN(1)+1,CLN(2)+1))
                  call ProjectFaceToNewPoints(SurfInfo(eIDRight) % facePatches(SideIDR), CLN(1), NodalStorage(CLN(1)) % xCGL, &
                                                                                         CLN(2), NodalStorage(CLN(2)) % xCGL, faceCL)
                  call SurfInfo(eIDRight) % facePatches(SideIDR) % Destruct()
                  call SurfInfo(eIDRight) % facePatches(SideIDR) % Construct(NodalStorage(CLN(1)) % xCGL,&
                                                                             NodalStorage(CLN(2)) % xCGL,faceCL)
                  deallocate(faceCL)
               end if

            case (HMESH_BOUNDARY)
               eIDLeft  = f % elementIDs(1)
               SideIDL  = f % elementSide(1)
               NSurfL   = SurfInfo(eIDLeft)  % facePatches(SideIDL) % noOfKnots - 1

               if     (SurfInfo(eIDLeft) % IsHex8 .or. all(NSurfL == 1)) cycle

               if (self % anisotropic  .and. (.not. self % meshIs2D) ) then
                  CLN = bfOrder(f % zone)
               else
                  CLN(1) = f % NfLeft(1)
                  CLN(2) = f % NfLeft(2)
               end if
!
!              Adapt the curved face order to the polynomial order
!              ---------------------------------------------------
               if ( any(CLN < NSurfL) ) then
                  allocate(faceCL(1:3,CLN(1)+1,CLN(2)+1))
                  call ProjectFaceToNewPoints(SurfInfo(eIDLeft) % facePatches(SideIDL), CLN(1), NodalStorage(CLN(1)) % xCGL, &
                                                                                        CLN(2), NodalStorage(CLN(2)) % xCGL, faceCL)
                  call SurfInfo(eIDLeft) % facePatches(SideIDL) % Destruct()
                  call SurfInfo(eIDLeft) % facePatches(SideIDL) % Construct(NodalStorage(CLN(1)) % xCGL, &
                                                                            NodalStorage(CLN(2)) % xCGL,faceCL)
                  deallocate(faceCL)
               end if

            case (HMESH_MPI)
               eID = maxval(f % elementIDs)
               thisSide = maxloc(f % elementIDs, dim=1)
               side = f % elementSide(thisSide)
               NSurf = SurfInfo(eID) % facePatches(side) % noOfKnots - 1

               if ( SurfInfo(eID) % IsHex8 .or. all(NSurf == 1) ) cycle

               if (self % elements(eID) % faceSide(side) == LEFT) then
                  CLN(1) = f % NfLeft(1)  
                  CLN(2) = f % NfLeft(2)  
               else
                  CLN(1) = f % NfRight(1)  
                  CLN(2) = f % NfRight(2)  
               end if

               if ( side .eq. 2 ) then    ! Right faces need to be rotated
                  select case ( f % rotation )
                  case ( 1, 3, 4, 6 ) ! Local x and y axis are perpendicular 
                     if (CLN(1) /= CLN(2)) then
                        buffer = CLN(1)
                        CLN(1) = CLN(2)
                        CLN(2) = buffer
                     end if
                  end select
               end if

               if ( any(CLN < NSurf) ) then
                  allocate(faceCL(1:3,CLN(1)+1,CLN(2)+1))
                  call ProjectFaceToNewPoints(SurfInfo(eID) % facePatches(side), CLN(1), NodalStorage(CLN(1)) % xCGL, &
                                                                                 CLN(2), NodalStorage(CLN(2)) % xCGL, faceCL)
                  call SurfInfo(eID) % facePatches(side) % Destruct()
                  call SurfInfo(eID) % facePatches(side) % Construct(NodalStorage(CLN(1)) % xCGL, &
                                                                     NodalStorage(CLN(2)) % xCGL,faceCL)
                  deallocate(faceCL)
               end if
            end select
            end associate
         end do
         safedeallocate (bfOrder)

!
!        ----------------------------
!        Construct elements' geometry
!        ----------------------------
!
         allocate(hex8Map)
         call hex8Map % constructWithCorners(corners)
         allocate(genHexMap)

         do ii=1, num_of_elements
            if ( present(elementList) ) then
               eID = elementList(ii)
            else
               eID = ii
            end if

            if (SurfInfo(eID) % IsHex8) then
               call hex8Map % setCorners(SurfInfo(eID) % corners)
               hexMap => hex8Map
            else
               CALL genHexMap % destruct()
               CALL genHexMap % constructWithFaces(SurfInfo(eID) % facePatches)

               hexMap => genHexMap
            end if

            call self % elements(eID) % ConstructGeometry (hexMap)

         end do
!
!        -------------------------
!        Construct faces' geometry
!        -------------------------
!
         do ii=1, num_of_faces
            if ( present(facesList) ) then
               fID = facesList(ii)
            else
               fID = ii
            end if

            associate(f => self % faces(fID))
            select case(f % faceType)
            case(HMESH_INTERIOR, HMESH_BOUNDARY)
               associate(eL => self % elements(f % elementIDs(1)))
               call f % geom % construct(f % Nf, f % NelLeft, f % NfLeft, eL % Nxyz, &
                                         NodalStorage(f % Nf), NodalStorage(eL % Nxyz), &
                                         eL % geom, eL % hexMap, f % elementSide(1), &
                                         f % projectionType(1), 1, 0 )
               end associate

            case(HMESH_MPI)
               side = maxloc(f % elementIDs, dim=1)

               select case (side)
               case(1)
                  Nel = f % NelLeft
                  Nelf = f % NfLeft
                  rot = 0

               case(2)
                  Nel = f % NelRight
                  Nelf = f % NfRight
                  rot = f % rotation

               end select

               associate(e => self % elements(f % elementIDs(side)))
               call f % geom % construct(f % Nf, Nelf, Nel, e % Nxyz, &
                                         NodalStorage(f % Nf), NodalStorage(e % Nxyz), &
                                         e % geom, e % hexMap, f % elementSide(side), &
                                         f % projectionType(side), side, rot)

               end associate

            end select
            end associate

         end do
!
!        ------------------------------------------------------
!        Compute the faces minimum orthogonal distance estimate
!        ------------------------------------------------------
!
         do ii = 1, num_of_faces
            if ( present(facesList) ) then
               fID = facesList(ii)
            else
               fID = ii
            end if
            associate(f => self % faces(fID))
            select case(f % faceType)
            case(HMESH_INTERIOR)
               f % geom % h = min(minval(self % elements(f % elementIDs(1)) % geom % jacobian), &
                                  minval(self % elements(f % elementIDs(2)) % geom % jacobian)) &
                        / maxval(f % geom % jacobian)
            case(HMESH_BOUNDARY)
               f % geom % h = minval(self % elements(f % elementIDs(1)) % geom % jacobian) &
                        / maxval(f % geom % jacobian)
            case(HMESH_MPI)
               f % geom % h = minval(self % elements(maxval(f % elementIDs)) % geom % jacobian) &
                        / maxval(f % geom % jacobian)
            end select
            end associate
         end do

         if ( MPI_Process % doMPIAction ) then
            call CommunicateMPIFaceMinimumDistance(self)
         end if

!
!        ---------
!        Finish up
!        ---------
!
         deallocate (SurfInfo)
         CALL hex8Map % destruct()
         DEALLOCATE(hex8Map)
         CALL genHexMap % destruct()
         DEALLOCATE(genHexMap)

      end subroutine HexMesh_ConstructGeometry

      subroutine CommunicateMPIFaceMinimumDistance(self)
         implicit none
         class(HexMesh) :: self
!
!        ---------------
!        Local variables
!        ---------------
!
#ifdef _HAS_MPI_
         integer  :: no_of_max_faces, i, fID, mpifID, ierr
         real(kind=RP), allocatable  :: hsend(:,:)
         real(kind=RP), allocatable  :: hrecv(:,:)
         integer  :: sendReq(MPI_Process % nProcs)
         integer  :: recvReq(MPI_Process % nProcs)
         integer  :: sendSt(MPI_STATUS_SIZE, MPI_Process % nProcs)
         integer  :: recvSt(MPI_STATUS_SIZE, MPI_Process % nProcs)

         if ( .not. self % MPIfaces % Constructed ) return
!
!        Get the maximum number of faces
!        -------------------------------
         no_of_max_faces = 0
         do i = 1, MPI_Process % nProcs
            no_of_max_faces = max(no_of_max_faces, self % MPIfaces % faces(i) % no_of_faces)
         end do
!
!        Allocate an array to store the distances
!        ----------------------------------------
         allocate(hsend(no_of_max_faces, MPI_Process % nProcs))
         allocate(hrecv(no_of_max_faces, MPI_Process % nProcs))
!
!        Perform the receive calls
!        -------------------------
         do i = 1, MPI_Process % nProcs
            if ( self % MPIfaces % faces(i) % no_of_faces .le. 0 ) then
               recvReq(i) = MPI_REQUEST_NULL
               cycle

            end if

            call mpi_irecv(hrecv(:,i), self % MPIfaces % faces(i) % no_of_faces, MPI_DOUBLE, i-1, MPI_ANY_TAG, &
                           MPI_COMM_WORLD, recvReq(i), ierr)

         end do

!
!        Gather the distances to send
!        ------------------------------
         do i = 1, MPI_Process % nProcs
            if ( self % MPIfaces % faces(i) % no_of_faces .le. 0 ) cycle

            do mpifID = 1, self % MPIfaces % faces(i) % no_of_faces
               fID = self % MPIfaces % faces(i) % faceIDs(mpifID)
               hsend(mpifID,i) = self % faces(fID) % geom % h
            end do
         end do
!
!        Send the distances
!        ------------------
         do i = 1, MPI_Process % nProcs
            if ( self % MPIfaces % faces(i) % no_of_faces .le. 0 ) then
               sendReq(i) = MPI_REQUEST_NULL
               cycle

            end if

            call mpi_isend(hsend(:,i), self % MPIfaces % faces(i) % no_of_faces, MPI_DOUBLE, i-1, DEFAULT_TAG, &
                           MPI_COMM_WORLD, sendReq(i), ierr)

         end do

         call mpi_waitall(MPI_Process % nProcs, recvReq, recvSt, ierr)
!
!        Collect the distances
!        ---------------------
         do i = 1, MPI_Process % nProcs
            if ( self % MPIfaces % faces(i) % no_of_faces .le. 0 ) cycle

            do mpifID = 1, self % MPIfaces % faces(i) % no_of_faces
               fID = self % MPIfaces % faces(i) % faceIDs(mpifID)
               self % faces(fID) % geom % h = min(self % faces(fID) % geom % h, hrecv(mpifID,i))
            end do
         end do

         call mpi_waitall(MPI_Process % nProcs, sendReq, sendSt, ierr)

         deallocate(hsend, hrecv)

#endif
      end subroutine CommunicateMPIFaceMinimumDistance

      subroutine HexMesh_DefineAsBoundaryFaces(self)
         implicit none
         class(HexMesh) :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: fID

         do fID = 1, size(self % faces)
            if ( self % faces(fID) % faceType .eq. HMESH_UNDEFINED ) then
               if ( trim(self % faces(fID) % boundaryName) .ne. emptyBCName ) then
                  self % faces(fID) % faceType = HMESH_BOUNDARY
               end if
            end if
         end do

      end subroutine HexMesh_DefineAsBoundaryFaces
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------
!     Export mesh to a hmesh file
!     ---------------------------
      subroutine HexMesh_Export(self, fileName)
         use SolutionFile
         use MPI_Process_Info
         implicit none
         class(HexMesh),   intent(in)     :: self
         character(len=*), intent(in)     :: fileName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: fid, eID
         integer(kind=AddrInt)         :: pos
         character(len=LINE_LENGTH)    :: meshName
         real(kind=RP), parameter      :: refs(NO_OF_SAVED_REFS) = 0.0_RP

!
!        Create file: it will be contained in ./MESH
!        -------------------------------------------
         meshName = "./MESH/" // trim(removePath(getFileName(fileName))) // ".hmesh"
         call CreateNewSolutionFile( trim(meshName), MESH_FILE, self % nodeType, &
                                     self % no_of_allElements, 0, 0.0_RP, refs)
!
!        Introduce all element nodal coordinates
!        ---------------------------------------
         fID = putSolutionFileInWriteDataMode(trim(meshName))
         do eID = 1, self % no_of_elements
            associate(e => self % elements(eID))
            pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + 3_AddrInt*e % offsetIO*SIZEOF_RP
            call writeArray(fID, e % geom % x(:,0:e%Nxyz(1),0:e%Nxyz(2),0:e%Nxyz(3))*Lref, position=pos)
            end associate
         end do
         close(fid)
!
!        Close the file
!        --------------
         call SealSolutionFile(trim(meshName))

      end subroutine HexMesh_Export
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_ExportBoundaryMesh(self,fileName)
         implicit none
         !-arguments--------------------------------
         class(HexMesh),   intent(in)     :: self
         character(len=*), intent(in)     :: fileName          !<  Name of file containing polynomial orders to initialize
         !-local-variables--------------------------
         integer                          :: fd       ! File unit
         integer                          :: zoneID, zfID, fID
         character(len=LINE_LENGTH)       :: bMeshName
         !------------------------------------------


!
!        Create file: it will be contained in ./MESH
!        -------------------------------------------
         bMeshName = "./MESH/" // trim(removePath(getFileName(fileName))) // ".bmesh"

         open( newunit = fd , file = trim(bMeshName), action = 'write')

!
!        Write file
!        ----------
         write(fd,*) size (self % zones)

         do zoneID = 1, size (self % zones)

            write (fd,*) self % zones(zoneID) % Name
            write (fd,*) self % zones(zoneID) % no_of_faces

            do zfID=1, self % zones(zoneID) % no_of_faces
               fID = self % zones(zoneID) % faces(zFID)
               write(fd,'(I8)',advance='no') self % faces(fID) % elementIDs(1)
            end do
            write(fd,*)

            do zfID=1, self % zones(zoneID) % no_of_faces
               fID = self % zones(zoneID) % faces(zFID)
               write(fd,'(I8)',advance='no') self % faces(fID) % elementSide(1)
            end do
            write(fd,*)

         end do

         close (fd)
      end subroutine HexMesh_ExportBoundaryMesh

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------------
!     Export mesh orders to a file
!     ----------------------------
      subroutine HexMesh_ExportOrders(self,fileName)
         implicit none
         !------------------------------------------
         class(HexMesh),   intent(in)     :: self
         character(len=*), intent(in)     :: fileName          !<  Name of file containing polynomial orders to initialize
         !------------------------------------------
         integer                          :: fd       ! File unit
         integer                          :: k
         character(len=LINE_LENGTH)       :: OrderFileName
         !------------------------------------------

!
!        Create file: it will be contained in ./MESH
!        -------------------------------------------
         OrderFileName = "./MESH/" // trim(removePath(getFileName(fileName))) // ".omesh"


         open( newunit = fd , FILE = TRIM(OrderFileName), ACTION = 'write')

            write(fd,*) size(self % elements)

            do k=1, size(self % elements)
               write(fd,*) self % elements(k) % Nxyz
            end do

         close (fd)

      end subroutine HexMesh_ExportOrders
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!     ************************************************************************
!           Save solution subroutine for the Navier-Stokes solver. It saves
!        the state vector (Q), and optionally the gradients.
!     ************************************************************************
!
     subroutine HexMesh_SaveSolution(self, iter, time, name, saveGradients)
         use SolutionFile
         use MPI_Process_Info
         implicit none
         class(HexMesh)                         :: self
         integer,             intent(in)        :: iter
         real(kind=RP),       intent(in)        :: time
         character(len=*),    intent(in)        :: name
         logical,             intent(in)        :: saveGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                          :: fid, eID, padding
         integer(kind=AddrInt)            :: pos
         real(kind=RP)                    :: refs(NO_OF_SAVED_REFS)
         real(kind=RP), allocatable       :: Q(:,:,:,:)
#if (!defined(NAVIERSTOKES))
         logical                          :: computeGradients = .true.
#endif
!
!        Gather reference quantities
!        ---------------------------
#if defined(NAVIERSTOKES)
         refs(GAMMA_REF) = thermodynamics % gamma
         refs(RGAS_REF)  = thermodynamics % R
         refs(RHO_REF)   = refValues      % rho
         refs(V_REF)     = refValues      % V
         refs(T_REF)     = refValues      % T
         refs(MACH_REF)  = dimensionless  % Mach
#elif defined(INCNS)
         refs(GAMMA_REF) = 0.0_RP
         refs(RGAS_REF)  = 0.0_RP
         refs(RHO_REF)   = refValues      % rho
         refs(V_REF)     = refValues      % V
         refs(T_REF)     = 0.0_RP
         refs(MACH_REF)  = 0.0_RP
#else
         refs = 0.0_RP
#endif
!
!        Create new file
!        ---------------
         if ( saveGradients .and. computeGradients) then
            call CreateNewSolutionFile(trim(name),SOLUTION_AND_GRADIENTS_FILE, &
                                       self % nodeType, self % no_of_allElements, iter, time, refs)
            padding = NCONS + 3*NGRAD
         else
            call CreateNewSolutionFile(trim(name),SOLUTION_FILE, self % nodeType, &
                                       self % no_of_allElements, iter, time, refs)
            padding = NCONS
         end if
!
!        Write arrays
!        ------------
         fID = putSolutionFileInWriteDataMode(trim(name))
         do eID = 1, self % no_of_elements
            associate( e => self % elements(eID) )

            allocate(Q(NCONS, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
#ifdef FLOW
            Q(1:NCONS,:,:,:)  = e % storage % Q
#ifdef MULTIPHASE
            Q(IMP,:,:,:) = e % storage % Q(IMP,:,:,:) + e % storage % Q(IMC,:,:,:)*e % storage % mu(1,:,:,:)
#endif
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
            Q(NCONS,:,:,:) = e % storage % c(1,:,:,:)
#endif

            pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + padding*e % offsetIO * SIZEOF_RP
            call writeArray(fid, Q, position=pos)

            deallocate(Q)
            if ( saveGradients .and. computeGradients ) then

               allocate(Q(NGRAD,0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))

#ifdef FLOW
               Q(1:NCONS,:,:,:) = e % storage % U_x
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
               Q(NGRAD,:,:,:) = e % storage % c_x(1,:,:,:)
#endif
               write(fid) Q

#ifdef FLOW
               Q(1:NCONS,:,:,:) = e % storage % U_y
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
               Q(NGRAD,:,:,:) = e % storage % c_y(1,:,:,:)
#endif
               write(fid) Q

#ifdef FLOW
               Q(1:NCONS,:,:,:) = e % storage % U_z
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
               Q(NGRAD,:,:,:) = e % storage % c_z(1,:,:,:)
#endif
               write(fid) Q

               deallocate(Q)
            end if

            end associate
         end do
         close(fid)
!
!        Close the file
!        --------------
         call SealSolutionFile(trim(name))

      end subroutine HexMesh_SaveSolution

#if defined(NAVIERSTOKES)
      subroutine HexMesh_SaveStatistics(self, iter, time, name, saveGradients)
         use SolutionFile
         implicit none
         class(HexMesh),      intent(in)        :: self
         integer,             intent(in)        :: iter
         real(kind=RP),       intent(in)        :: time
         character(len=*),    intent(in)        :: name
         logical,             intent(in)        :: saveGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                          :: fid, eID
         integer                          :: no_stat_s
         integer(kind=AddrInt)            :: pos
         real(kind=RP)                    :: refs(NO_OF_SAVED_REFS) 
         real(kind=RP), allocatable       :: Q(:,:,:,:)
!
!        Gather reference quantities
!        ---------------------------
         refs(GAMMA_REF) = thermodynamics % gamma
         refs(RGAS_REF)  = thermodynamics % R
         refs(RHO_REF)   = refValues      % rho
         refs(V_REF)     = refValues      % V
         refs(T_REF)     = refValues      % T
         refs(MACH_REF)  = dimensionless  % Mach
         refs(RE_REF)    = dimensionless  % Re

!        Create new file
!        ---------------
         call CreateNewSolutionFile(trim(name),STATS_FILE, self % nodeType, self % no_of_allElements, iter, time, refs)
!
!        Write arrays
!        ------------
         fID = putSolutionFileInWriteDataMode(trim(name))
         do eID = 1, self % no_of_elements
            associate( e => self % elements(eID) )
            pos = POS_INIT_DATA + (e % globID-1)*5_AddrInt*SIZEOF_INT + no_of_stats_variables*e % offsetIO*SIZEOF_RP
            no_stat_s = 9
            call writeArray(fid, e % storage % stats % data(1:no_stat_s,:,:,:), position=pos)
            allocate(Q(NCONS, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
            ! write(fid) e%storage%stats%data(7:,:,:,:)
            Q(1:NCONS,:,:,:) = e % storage % stats % data(no_stat_s+1:no_stat_s+NCONS,:,:,:)
            write(fid) Q
            deallocate(Q)
            if ( saveGradients .and. computeGradients ) then
               allocate(Q(NGRAD,0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
               ! UX
               Q(1:NGRAD,:,:,:) = e % storage % stats % data(no_stat_s+NCONS+1:no_stat_s+NCONS+NGRAD,:,:,:)
               write(fid) Q
               ! UY
               Q(1:NGRAD,:,:,:) = e % storage % stats % data(no_stat_s+NCONS+1+NGRAD:no_stat_s+NCONS+2*NGRAD,:,:,:)
               write(fid) Q
               ! UZ
               Q(1:NGRAD,:,:,:) = e % storage % stats % data(no_stat_s+NCONS+1+2*NGRAD:,:,:,:)
               write(fid) Q
               deallocate(Q)
            end if
            end associate
         end do
         close(fid)
!
!        Close the file
!        --------------
         call SealSolutionFile(trim(name))

      end subroutine HexMesh_SaveStatistics

      subroutine HexMesh_ResetStatistics(self)
         implicit none
         class(HexMesh)       :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID

         do eID = 1, self % no_of_elements
            self % elements(eID) % storage % stats % data = 0.0_RP
         end do

      end subroutine HexMesh_ResetStatistics
#endif
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------------------------
!     Subroutine to load a solution for restart using the information in the control file
!     -----------------------------------------------------------------------------------
      subroutine HexMesh_LoadSolutionForRestart( self, controlVariables, initial_iteration, initial_time, loadFromNSSA ) 
         use mainKeywordsModule, only: restartFileNameKey
         use FileReaders       , only: ReadOrderFile
         implicit none
         !-arguments-----------------------------------------------
         class(HexMesh)                       :: self
         type(FTValueDictionary), intent(in)  :: controlVariables
         logical                , intent(in)  :: loadFromNSSA
         integer                , intent(out) :: initial_iteration
         real(kind=RP)          , intent(out) :: initial_time
         !-local-variables-----------------------------------------
         character(len=LINE_LENGTH)           :: fileName
         character(len=LINE_LENGTH)           :: orderFileName
         integer, allocatable                 :: Nx(:), Ny(:), Nz(:)
         type(HexMesh), target                :: auxMesh
         integer                              :: NDOF, eID
         logical                              :: with_gradients
#if (!defined(NAVIERSTOKES)) || (!defined(INCNS)) || (!defined(MULTIPHASE))
         logical                          :: computeGradients = .true.
#endif
         !---------------------------------------------------------

         fileName = controlVariables % stringValueForKey(restartFileNameKey,requestedLength = LINE_LENGTH)

!
!        *****************************************************
!        The restart polynomial orders are different to self's
!        *****************************************************
!
         if ( controlVariables % containsKey("restart polorder")        .or. &
              controlVariables % containsKey("restart polorder file") ) then


!           Read the polynomial order of the solution to be loaded
!           ------------------------------------------------------

            if ( controlVariables % containsKey("restart polorder") ) then
               allocate ( Nx(self % no_of_allElements) , Ny(self % no_of_allElements) , Nz(self % no_of_allElements) )
               Nx = controlVariables % integerValueForKey ("restart polorder")
               Ny = Nx
               Nz = Nx
            elseif ( controlVariables % containsKey("restart polorder file") ) then
               orderFileName = controlVariables % stringValueForKey("restart polorder file",requestedLength = LINE_LENGTH)
               call ReadOrderFile(trim(orderFileName), Nx, Ny, Nz)
            end if

!           Construct an auxiliar mesh to read the solution
!           -----------------------------------------------

            auxMesh % nodeType = self % nodeType
            auxMesh % no_of_elements = self % no_of_elements
            auxMesh % no_of_allElements = self % no_of_allElements
            allocate ( auxMesh % elements (self % no_of_elements) )
            allocate ( auxMesh % Nx (self % no_of_elements) )
            allocate ( auxMesh % Ny (self % no_of_elements) )
            allocate ( auxMesh % Nz (self % no_of_elements) )

            NDOF = 0
            do eID = 1, self % no_of_elements
               associate ( e_aux => auxMesh % elements(eID), &
                           e     =>    self % elements(eID) )
               auxMesh % Nx(eID) = Nx (e % globID)
               auxMesh % Ny(eID) = Ny (e % globID)
               auxMesh % Nz(eID) = Nz (e % globID)
               e_aux % globID = e % globID
               e_aux % Nxyz = [Nx(e % globID) , Ny(e % globID) , Nz(e % globID)]
               NDOF = NDOF + (Nx(e % globID) + 1) * (Ny(e % globID) + 1) * (Nz(e % globID) + 1)             
               end associate
            end do

            call auxMesh % PrepareForIO
            call auxMesh % AllocateStorage (NDOF, controlVariables,computeGradients,.FALSE.)
            call auxMesh % storage % pointStorage
            do eID = 1, auxMesh % no_of_elements
               auxMesh % elements(eID) % storage => auxMesh % storage % elements(eID)
            end do

!           Read the solution in the auxiliar mesh and interpolate to current mesh
!           ----------------------------------------------------------------------
            
            call auxMesh % LoadSolution ( fileName, initial_iteration, initial_time , with_gradients, loadFromNSSA=loadFromNSSA)
            do eID=1, self % no_of_elements
               call auxMesh % storage % elements (eID) % InterpolateSolution (self % storage % elements(eID), auxMesh % nodeType , with_gradients)
            end do

!           Clean up
!           --------

            do eID = 1, auxMesh % no_of_elements
               call auxMesh % elements(eID) % storage % destruct
            end do
            call auxMesh % storage % destruct
            deallocate (auxMesh % elements)

            if ( controlVariables % containsKey("get discretization error of") ) call GetDiscretizationError(self,controlVariables)
!
!        *****************************************************
!        The restart polynomial orders are the same as in self
!        *****************************************************
!
         else
            call self % LoadSolution ( fileName, initial_iteration, initial_time, loadFromNSSA=loadFromNSSA )
         end if

      end subroutine HexMesh_LoadSolutionForRestart
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------
!     Subroutine to load a solution (*.hsol) in self
!     ----------------------------------------------
      subroutine HexMesh_LoadSolution( self, fileName, initial_iteration, initial_time , with_gradients, loadFromNSSA)
         IMPLICIT NONE
         CLASS(HexMesh)                  :: self
         character(len=*)                :: fileName
         integer           , intent(out) :: initial_iteration
         real(kind=RP)     , intent(out) :: initial_time
         logical, optional , intent(out) :: with_gradients
         logical, optional , intent(in)  :: loadFromNSSA
         
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                        :: fID, eID, fileType, no_of_elements, flag, nodetype
         integer                        :: padding, pos
         integer                        :: Nxp1, Nyp1, Nzp1, no_of_eqs, array_rank
         real(kind=RP), allocatable     :: Q(:,:,:,:)
         character(len=SOLFILE_STR_LEN) :: rstName
         logical          :: gradients
         logical                        :: NS_from_NSSA
         
         gradients = .FALSE.
         if (present(loadFromNSSA)) then
             NS_from_NSSA = loadFromNSSA
         else
             NS_from_NSSA = .FALSE.
         end if 
!
!        Get the file title
!        ------------------
         rstName = getSolutionFileName(trim(fileName))
!
!        Get the file type
!        -----------------
         fileType = getSolutionFileType(trim(fileName))

         select case (fileType)
         case(MESH_FILE)
            print*, "The selected restart file is a mesh file"
            errorMessage(STD_OUT)
            stop

         case(SOLUTION_FILE)
            padding = 1*NCONS

         case(SOLUTION_AND_GRADIENTS_FILE)
#ifdef SPALARTALMARAS
            padding = NCONS + 3 * NGRAD
            gradients = .TRUE.
#else
            padding = NCONS + 3 * NGRAD
            gradients = .TRUE.

#endif
         case(STATS_FILE)
            print*, "The selected restart file is a statistics file"
            errorMessage(STD_OUT)
            stop
         case default
            print*, "Unknown restart file format"
            errorMessage(STD_OUT)
            stop
         end select
         if (NS_from_NSSA) then
             if (gradients) then
                 ! add 1 as NNSA has one more NCONS, and 3 as has one more NGRAD, the whole will be padding = (NCONS + 1) + 3*(NGRAD+1)
                 padding = padding + 1 + 3
             else
                 padding = padding + 1
             end if
         end if
!
!        Get the node type
!        -----------------
         nodeType = getSolutionFileNodeType(trim(fileName))

         if ( nodeType .ne. self % nodeType ) then
            print*, "Solution file uses a different discretization nodes than the mesh."
            errorMessage(STD_OUT)
         end if
!
!        Read the number of elements
!        ---------------------------
         no_of_elements = getSolutionFileNoOfElements(trim(fileName))

         if ( no_of_elements .ne. self % no_of_allElements ) then
            write(STD_OUT,'(A,A)') "The number of elements stored in the restart file ", &
                                   "do not match that of the mesh file"
            errorMessage(STD_OUT)
            stop
         end if
!
!        Read the initial iteration and time
!        -----------------------------------
         call getSolutionFileTimeAndITeration(trim(fileName), initial_iteration, initial_time)
!
!        Read the terminator indicator
!        -----------------------------
         flag = getSolutionFileDataInitFlag(trim(fileName))

         if ( flag .ne. BEGINNING_DATA ) then
            print*, "Beginning data flag was not found in the file."
            errorMessage(STD_OUT)
            stop
         end if
!
!        Read elements data
!        ------------------
         fID = putSolutionFileInReadDataMode(trim(fileName))
         do eID = 1, size(self % elements)
            associate( e => self % elements(eID) )
            pos = POS_INIT_DATA + (e % globID-1)*5*SIZEOF_INT + padding*e % offsetIO*SIZEOF_RP
            read(fID, pos=pos) array_rank
            read(fID) no_of_eqs, Nxp1, Nyp1, Nzp1
            if (      ((Nxp1-1) .ne. e % Nxyz(1)) &
                 .or. ((Nyp1-1) .ne. e % Nxyz(2)) &
                 .or. ((Nzp1-1) .ne. e % Nxyz(3)) &
                 .or. (no_of_eqs .ne. NCONS )       ) then
               write(STD_OUT,'(A,I0,A)') "Error reading restart file: wrong dimension for element "&
                                           ,eID,"."

               write(STD_OUT,'(A,I0,A,I0,A,I0,A)') "Element dimensions: ", e % Nxyz(1), &
                                                                     " ,", e % Nxyz(2), &
                                                                     " ,", e % Nxyz(3), &
                                                                     "."

               write(STD_OUT,'(A,I0,A,I0,A,I0,A)') "Restart dimensions: ", Nxp1-1, &
                                                                     " ,", Nyp1-1, &
                                                                     " ,", Nzp1-1, &
                                                                     "."

               errorMessage(STD_OUT)
               stop
            end if

            allocate(Q(NCONS, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
            read(fID) Q
#ifdef FLOW
            e % storage % QNS = Q(1:NCONS,:,:,:)
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
            e % storage % c(1,:,:,:) = Q(NCONS,:,:,:)
#endif

            deallocate(Q)

            if (gradients) then
               allocate(Q(NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
               read(fID) Q
#ifdef FLOW
               e % storage % U_x = Q(1:NCONS,:,:,:)
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
               e % storage % c_x(1,:,:,:) = Q(NGRAD,:,:,:)
#endif

               read(fID) Q
#ifdef FLOW
               e % storage % U_y = Q(1:NCONS,:,:,:)
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
               e % storage % c_y(1,:,:,:) = Q(NGRAD,:,:,:)
#endif
               read(fID) Q
#ifdef FLOW
               e % storage % U_z = Q(1:NCONS,:,:,:)
#endif
#if (defined(CAHNHILLIARD) && (!defined(FLOW)))
               e % storage % c_z(1,:,:,:) = Q(NGRAD,:,:,:)
#endif
               deallocate(Q)
            end if

           end associate
         end do
!
!        Close the file
!        --------------
         close(fID)

         if (present(with_gradients) ) with_gradients = gradients

      END SUBROUTINE HexMesh_LoadSolution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine GetDiscretizationError(mesh,controlVariables)
      implicit none
      !----------------------------------------------
      type(HexMesh) :: mesh
      type(FTValueDictionary) :: controlVariables
      !----------------------------------------------
      character(len=LINE_LENGTH) :: fileName
      type(HexMesh)              :: refMesh
      integer       :: NDOF, eID
      integer       :: initial_iteration
      real(kind=RP) :: initial_time
      !----------------------------------------------

      fileName = controlVariables % stringValueForKey("get discretization error of",LINE_LENGTH)

!     Construct an auxiliar mesh to read the solution
!     -----------------------------------------------

      refMesh % nodeType = mesh % nodeType
      refMesh % no_of_elements = mesh % no_of_elements
      allocate ( refMesh % elements (mesh % no_of_elements) )

      NDOF = 0
      do eID = 1, mesh % no_of_elements
         associate ( e_aux => refMesh % elements(eID), &
                     e     =>    mesh % elements(eID) )
         e_aux % globID = e % globID
         e_aux % Nxyz = e % Nxyz
         NDOF = NDOF + product(e % Nxyz + 1)
         end associate
      end do

      call refMesh % PrepareForIO
      call refMesh % AllocateStorage (NDOF, controlVariables,.FALSE.)

!     Read the solution in the auxiliar mesh and interpolate to current mesh
!     ----------------------------------------------------------------------

      call refMesh % LoadSolution ( fileName, initial_iteration, initial_time )

      refMesh % storage % Q = mesh % storage % Q - refMesh % storage % Q

      print*, '|disc_error|∞ = ', maxval(abs(refMesh % storage % Q))
      print*, '|disc_error|² = ', norm2(refMesh % storage % Q)

      call refMesh % SaveSolution(0, 0._RP, 'RESULTS/DiscError.hsol', .FALSE.)

!           Clean up
!           --------

      do eID = 1, refMesh % no_of_elements
         call refMesh % elements(eID) % storage % destruct
      end do
      call refMesh % storage % destruct
      deallocate (refMesh % elements)

   end subroutine GetDiscretizationError
!
!////////////////////////////////////////////////////////////////////////
!
!        AUXILIAR SUBROUTINES
!        --------------------
!
!////////////////////////////////////////////////////////////////////////
!
      logical function HexMesh_FindPointWithCoords(self, x, eID, xi, optionalElements)
         implicit none
         class(HexMesh), intent(in)         :: self
         real(kind=RP),    intent(in)       :: x(NDIM)
         integer,          intent(out)      :: eID
         real(kind=RP),    intent(out)      :: xi(NDIM)
         integer, optional,intent(in)       :: optionalElements(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: op_eID
         integer     :: zoneID, fID
         logical     :: success

         HexMesh_FindPointWithCoords = .false.
!
!        Search in optionalElements (if present)
!        ---------------------------------------
         if ( present(optionalElements) ) then
            do op_eID = 1, size(optionalElements)
               if ( optionalElements(op_eID) .eq. -1 ) cycle
               associate(e => self % elements(optionalElements(op_eID)))
               success = e % FindPointWithCoords(x, self % dir2D_ctrl, xi)
               if ( success ) then
                  eID = optionalElements(op_eID)
                  HexMesh_FindPointWithCoords = .true.
                  return
               end if
               end associate
            end do
         end if
!
!        Search in linear (not curved) mesh (faster and safer)
!        -----------------------------------------------------
         do eID = 1, self % no_of_elements
            success = self % elements(eID) % FindPointInLinElement(x, self % nodes)
            if ( success ) exit
         end do
!
!        If found in linear mesh, use FindPointWithCoords in that element and, if necessary, in neighbors...
!        ---------------------------------------------------------------------------------------------------
         if (eID <= self % no_of_elements) then
            success = self % FindPointWithCoordsInNeighbors(x, xi, eID, 2)
            if ( success ) then
               HexMesh_FindPointWithCoords = .true.
               return
            end if
         end if
!
!        As a last resource, search using FindPointWithCoords only in boundary elements
!        ------------------------------------------------------------------------------
         do zoneID=1, size(self % zones)
            do fID=1, self % zones(zoneID) % no_of_faces

               op_eID = self % faces ( self % zones(zoneID) % faces(fID) ) % elementIDs(1)
               success = self % elements (op_eID) % FindPointWithCoords(x, self % dir2D_ctrl, xi)
               if ( success ) then
                  HexMesh_FindPointWithCoords = .true.
                  return
               end if
            end do
         end do

      end function HexMesh_FindPointWithCoords
!
!////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------------------------------------------------------------------------------
!     HexMesh_FindPointWithCoordsInNeighbors:
!     This subroutine looks for a point (defined by coordinates) in the neighbor elements of a specific element
!     Note: For MPI, this routine ONLY checks in neighbors that are in the same partition...
!     ---------------------------------------------------------------------------------------------------------
      logical recursive function HexMesh_FindPointWithCoordsInNeighbors(self, x, xi, eID, depth)
         implicit none
         !-arguments--------------------------------------------------
         class(HexMesh), intent(in)     :: self
         real(kind=RP) , intent(in)     :: x(NDIM)
         real(kind=RP) , intent(out)    :: xi(NDIM)
         integer       , intent(inout)  :: eID
         integer       , intent(in)     :: depth
         !-local-variables--------------------------------------------
         logical :: success
         integer :: fID, nID, new_eID
         !------------------------------------------------------------

         success = self % elements(eID) % FindPointWithCoords(x, self % dir2D_ctrl, xi)
         if ( success ) then
            HexMesh_FindPointWithCoordsInNeighbors = .TRUE.
            return
         end if

         if (depth > 1) then
            do fID=1, FACES_PER_ELEMENT

               new_eID = mpi_partition % global2localeID (self % elements(eID) % Connection(fID) % globID)
               if (new_eID == 0) cycle
               success = self % FindPointWithCoordsInNeighbors(x, xi, new_eID, depth-1)
               if ( success ) then
                  HexMesh_FindPointWithCoordsInNeighbors = .TRUE.
                  return
               end if

            end do
         end if

      end function HexMesh_FindPointWithCoordsInNeighbors
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_ComputeWallDistances(self,facesList,elementList)
         implicit none
         class(HexMesh)     :: self
         integer, optional , intent(in)    :: facesList(:)
         integer, optional , intent(in)    :: elementList(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: eID, ii, i, j, k, no_of_wallDOFS, num_of_elems, num_of_faces
         real(kind=RP) :: currentDistance, minimumDistance
         integer       :: fID
         real(kind=RP) :: xP(NDIM)
         real(kind=RP), allocatable    :: Xwall(:,:)
!
!        Gather all walls coordinates
!        ----------------------------
         call HexMesh_GatherAllWallCoordinates(self, no_of_wallDOFS, Xwall)
!
!        Get the minimum distance to each elements nodal degree of freedom
!        -----------------------------------------------------------------
         if ( present(elementList) ) then
            num_of_elems = size (elementList)
         else
            num_of_elems = self % no_of_elements
         end if
         do ii = 1, num_of_elems
            if ( present(elementList) ) then
               eID = elementList (ii)
            else
               eID = ii
            end if
            associate(e => self % elements(eID))
            allocate(e % geom % dWall(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))

            do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               xP = e % geom % x(:,i,j,k)

               minimumDistance = HUGE(1.0_RP)
               do fID = 1, no_of_wallDOFS
                  currentDistance = sum(POW2(xP - Xwall(:,fID)))
                  minimumDistance = min(minimumDistance, currentDistance)
               end do

               e % geom % dWall(i,j,k) = sqrt(minimumDistance)

            end do                  ; end do                ; end do
            end associate
         end do
!
!        Get the minimum distance to each face nodal degree of freedom
!        -------------------------------------------------------------
         if ( present(facesList) ) then
            num_of_faces = size (facesList)
         else
            num_of_faces = size (self % faces)
         end if
         do ii = 1, num_of_faces
            if ( present(facesList) ) then
               eID = facesList (ii)
            else
               eID = ii
            end if

            associate(fe => self % faces(eID))
            allocate(fe % geom % dWall(0:fe % Nf(1), 0:fe % Nf(2)))

            do j = 0, fe % Nf(2) ; do i = 0, fe % Nf(1)
               xP = fe % geom % x(:,i,j)

               minimumDistance = HUGE(1.0_RP)
               do fID = 1, no_of_wallDOFS
                  currentDistance = sum(POW2(xP - Xwall(:,fID)))
                  minimumDistance = min(minimumDistance, currentDistance)
               end do

               fe % geom % dWall(i,j) = sqrt(minimumDistance)

            end do                ; end do
            end associate
         end do

         deallocate(Xwall)

      end subroutine HexMesh_ComputeWallDistances

      subroutine HexMesh_GatherAllWallCoordinates(self, no_of_wallDOFS, Xwall)
         implicit none
         class(HexMesh),              intent(in)  :: self
         integer,                     intent(out) :: no_of_wallDOFS
         real(kind=RP),  allocatable, intent(out) :: Xwall(:,:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: no_of_localWallDOFS
         integer                    :: no_of_wallDOFS_perProcess(MPI_Process % nProcs)
         integer                    :: displ(MPI_Process % nProcs)
         real(kind=RP), allocatable :: localXwall(:,:)
         integer                    :: zID, ierr, fID, zonefID, i, j, displacement

         no_of_localWallDOFS = 0
         do zID = 1, size(self % zones)
            if ( (trim(BCs(zID) % bc % bcType) .ne. "noslipwall") ) then
               cycle
            end if

            do zonefID = 1, self % zones(zID) % no_of_faces
               fID = self % zones(zID) % faces(zonefID)
               no_of_localWallDOFS = no_of_localWallDOFS + product(self % faces(fID) % Nf + 1)
            end do
         end do
         allocate( localXwall(NDIM, no_of_localWallDOFS) )
!
!        Loop all faces to gather the local wall coordinates
!        ---------------------------------------------------
         no_of_localWallDOFS = 0
         do zID = 1, size(self % zones)
            if ( (trim(BCs(zID) % bc % bcType) .ne. "noslipwall") ) then
               cycle
            end if

            do zonefID = 1, self % zones(zID) % no_of_faces
               fID = self % zones(zID) % faces(zonefID)
               associate( f => self % faces(fID) )
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  no_of_localWalLDOFS = no_of_localWallDOFS + 1
                  localXWall(1:NDIM,no_of_localWallDOFS) = f % geom % x(1:NDIM,i,j)
               end do               ; end do
               end associate
            end do
         end do
!
!        ******************************************************************
!        MPI: To communicate each process walls, we use the allgatherv
!             function
!        ******************************************************************
!
         if ( MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
!
!           Perform a reduction to how many walls are in each process
!           ---------------------------------------------------------
            call mpi_allgather(no_of_localWallDOFS, 1, MPI_INT, no_of_wallDOFS_perProcess, 1, MPI_INT, MPI_COMM_WORLD, ierr)
!
!           Compute the displacements
!           -------------------------
            displ(1) = 0
            do zID = 1, MPI_Process % nProcs - 1
               displ(zID+1) = displ(zID) + no_of_wallDOFS_perProcess(zID) * NDIM
            end do
!
!           Allocate the data
!           -----------------
            no_of_wallDOFS = sum(no_of_wallDOFS_perProcess)
            allocate( Xwall(NDIM, no_of_wallDOFS) )
!
!           Perform an allgatherv
!           ---------------------
            call mpi_allgatherv(localXwall, NDIM*no_of_localWallDOFS, MPI_DOUBLE, Xwall, NDIM*no_of_wallDOFS_perProcess, displ, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
#endif
         else
            no_of_wallDOFS = no_of_localWallDOFS
            allocate( Xwall(NDIM, no_of_wallDOFS) )
            Xwall = localXwall
         end if

         deallocate(localXwall)

      end subroutine HexMesh_GatherAllWallCoordinates
!
!////////////////////////////////////////////////////////////////////////
!
!        CONSTRUCT ZONES
!        ---------------
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_ConstructZones( self )
      implicit none
      class(HexMesh)          :: self

      call ConstructZones ( self % faces , self % zones )

      end subroutine HexMesh_ConstructZones

!
!///////////////////////////////////////////////////////////////////////
!
!  HexMesh_AllocateStorage:
!  Allocates the storage for the simulation
!  -> Storage specific to the analytical Jacobian is constructed by the corresponding class (AnalyticalJacobian.f90)
!
   subroutine HexMesh_AllocateStorage(self,NDOF,controlVariables,computeGradients,Face_Storage)
      implicit none
      !-----------------------------------------------------------
      class(HexMesh), target                 :: self
      integer                 , intent(in)   :: NDOF
      class(FTValueDictionary), intent(in)   :: controlVariables
      logical                 , intent(in)   :: computeGradients
      logical, optional       , intent(in)   :: Face_Storage
      !-----------------------------------------------------------
      integer :: bdf_order, eID, fID, RKSteps_num
      logical :: Face_St, FaceComputeQdot
      character(len=LINE_LENGTH) :: time_int
      character(len=LINE_LENGTH) :: mg_smoother
      !-----------------------------------------------------------

      if ( present(Face_Storage) ) then
         Face_St = Face_Storage
      else
         Face_St = .TRUE.
      end if
      FaceComputeQdot = controlVariables % containsKey("acoustic analogy")

      time_int = controlVariables % stringValueForKey("time integration",LINE_LENGTH)
      call toLower (time_int)

      if     ( controlVariables % containsKey("bdf order")) then
         bdf_order = controlVariables % integerValueForKey("bdf order")
         RKSteps_num = 0
      elseif ( trim(time_int) == "fas" ) then
         bdf_order = -1
         RKSteps_num = 0
        if ( controlVariables % containsKey("mg smoother")) then
          mg_smoother = controlVariables % stringValueForKey("mg smoother",LINE_LENGTH)
          call toLower (mg_smoother)
          if ( (trim(mg_smoother) .eq. "irk") .or. (trim(mg_smoother) .eq. "birk5") &
            .or. (trim(mg_smoother) .eq. "ilu") .or. (trim(mg_smoother) .eq. "sgs") ) then
            bdf_order = 1
            RKSteps_num = 0
          end if
        end if
      elseif ( trim(time_int) == "imex" ) then
         bdf_order = 1
#ifdef MULTIPHASE
         RKSteps_num = 3
!         RKSteps_num = 0
#else
         RKSteps_num = 0
#endif
      elseif ( trim(time_int) == "rosenbrock" ) then
         bdf_order = 0
         RKSteps_num = 0
      else
         bdf_order = -1
         RKSteps_num = 0
      end if

!     Construct global and elements' storage
!     --------------------------------------
      call self % storage % construct (NDOF, self % Nx, self % Ny, self % Nz, computeGradients, .FALSE., bdf_order, RKSteps_num )

!     Construct faces' storage
!     ------------------------
      if (Face_St) then
         do fID = 1, size(self % faces)
            associate ( f => self % faces(fID) )
            call f % storage(1) % Construct(NDIM, f % Nf, f % NelLeft , computeGradients, .FALSE., FaceComputeQdot)
            call f % storage(2) % Construct(NDIM, f % Nf, f % NelRight, computeGradients, .FALSE., FaceComputeQdot)
            end associate
         end do
      end if

!     Point element storage
!     ---------------------
      DO eID = 1, SIZE(self % elements)
         associate (e => self % elements(eID))
         e % hn = (e % geom % Volume / product(e % Nxyz + 1)) ** (1.0_RP / 3.0_RP)  ! Also compute h/p here
         e % storage => self % storage % elements(eID)
         end associate
      END DO

   end subroutine HexMesh_AllocateStorage

   subroutine HexMesh_SetStorageToEqn(self, which)
      implicit none
      class(HexMesh), target :: self
      integer, intent(in)    :: which
!
!     ---------------
!     Local variables
!     ---------------
!
      integer  :: off, ns, c, mu, nssa
      integer  :: eID, fID

      call GetStorageEquations(off, ns, c, mu, nssa)

      if ( which .eq. ns ) then
#ifdef FLOW
         self % storage % Q => self % storage % QNS
         self % storage % QDot => self % storage % QDotNS
         self % storage % PrevQ(1:,1:) => self % storage % PrevQNS(1:,1:)

         do eID = 1, self % no_of_elements
            call self % elements(eID) % storage % SetStorageToNS
         end do

         do fID = 1, size(self % faces)
            call self % faces(fID) % storage(1) % SetStorageToNS
            call self % faces(fID) % storage(2) % SetStorageToNS
         end do

      elseif ( which .eq. nssa ) then

         self % storage % Q => self % storage % QNS
         self % storage % QDot => self % storage % QDotNS
         self % storage % PrevQ(1:,1:) => self % storage % PrevQNS(1:,1:)

         do eID = 1, self % no_of_elements
            call self % elements(eID) % storage % SetStorageToNS
         end do

         do fID = 1, size(self % faces)
            call self % faces(fID) % storage(1) % SetStorageToNS
            call self % faces(fID) % storage(2) % SetStorageToNS
         end do


#endif


      elseif ( which .eq. c ) then
#if defined(CAHNHILLIARD)
         self % storage % Q => self % storage % c
         self % storage % QDot => self % storage % cDot
         self % storage % PrevQ => self % storage % PrevC

         do eID = 1, self % no_of_elements
            call self % elements(eID) % storage % SetStorageToCH_c
         end do

         do fID = 1, size(self % faces)
            call self % faces(fID) % storage(1) % SetStorageToCH_c
            call self % faces(fID) % storage(2) % SetStorageToCH_c
         end do
#endif
      elseif ( which .eq. mu ) then
#if defined(CAHNHILLIARD)
         self % storage % Q => self % storage % c
         self % storage % QDot => self % storage % cDot
         self % storage % PrevQ => self % storage % PrevC

         do eID = 1, self % no_of_elements
            call self % elements(eID) % storage % SetStorageToCH_mu
         end do

         do fID = 1, size(self % faces)
            call self % faces(fID) % storage(1) % SetStorageToCH_mu
            call self % faces(fID) % storage(2) % SetStorageToCH_mu
         end do
#endif
      end if

   end subroutine HexMesh_SetStorageToEqn
!
!///////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------------------------
!  Checks if the representation is conforming on a zone (Boundary)
!  ---------------------------------------------------------------
   function HexMesh_ConformingOnZone (self,zoneID) result(conforming)
      implicit none
      !-----------------------------------------------------------
      class(HexMesh), intent(in) :: self
      integer       , intent(in) :: zoneID
      logical                    :: conforming
      !-----------------------------------------------------------
      integer :: fIdx   ! Local face index in zone
      integer :: fID    ! Face index
      integer :: eID    ! Element index
      integer :: eSide  ! Side of the element in contact with boundary
      integer :: nFace  ! Counter for neighbor faces
      !-----------------------------------------------------------

      if (zoneID < lbound(self % zones,1) .or. zoneID > ubound(self % zones,1) ) ERROR stop 'HexMesh_ConformingOnZone :: Out of bounds'

      conforming = .TRUE.
      do fIdx = 1, self % zones(zoneID) % no_of_faces

         fID   = self % zones(zoneID) % faces(fIdx)
         eID   = self % faces(fID) % elementIDs(1)
         eSide = self % faces(fID) % elementSide(1)

         ! loop over the faces that are shared between boundary elements
         do nFace = 1, 4
            associate (f => self % faces ( self % elements(eID) % faceIDs (neighborFaces(nFace,eSide) ) ) )

            if (f % FaceType == HMESH_BOUNDARY) cycle

            if (any(f % NfLeft /= f % NfRight)) then
               conforming = .FALSE.
               return
            end if

            end associate
         end do

      end do

   end function
#if defined(INCNS) && defined(CAHNHILLIARD)
   subroutine HexMesh_ConvertDensityToPhaseField(self)
!
!     *************************************************************
!     Convert density to phase field only in element interior nodes
!     *************************************************************
!
      implicit none
      class(HexMesh),   intent(inout)  :: self
!
!     ---------------
!     Local variables
!     ---------------
!
      integer  :: eID, fID

      associate(rho1 => dimensionless % rho(1), &
                rho2 => dimensionless % rho(2))
      do eID = 1, self % no_of_elements
         associate(c => self % elements(eID) % storage % c, &
                   Q => self % elements(eID) % storage % QNS)
         c(1,:,:,:) = (-rho1 - rho2 + 2.0_RP * Q(INSRHO,:,:,:))/(rho2-rho1)
         end associate
      end do

      end associate
   end subroutine HexMesh_ConvertDensityToPhaseField

   subroutine HexMesh_ConvertPhaseFieldToDensity(self)
!
!     *************************************************************
!     Convert density to phase field only in element interior nodes
!     *************************************************************
!
      implicit none
      class(HexMesh),   intent(inout)  :: self
!
!     ---------------
!     Local variables
!     ---------------
!
      integer  :: eID, fID

      associate(rho1 => dimensionless % rho(1), &
                rho2 => dimensionless % rho(2))
      do eID = 1, self % no_of_elements
         associate(c => self % elements(eID) % storage % c, &
                   Q => self % elements(eID) % storage % QNS)
         Q(INSRHO,:,:,:) = 0.5_RP*(rho1*(1.0_RP-c(1,:,:,:)) + rho2*(1.0_RP + c(1,:,:,:)))
         end associate
      end do

      end associate

   end subroutine HexMesh_ConvertPhaseFieldToDensity
#endif
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------------
!  Adapts a mesh to new polynomial orders NNew
!  --------------------------------------------------------
   subroutine HexMesh_pAdapt (self, NNew, controlVariables)
      implicit none
      !-arguments-----------------------------------------
      class(HexMesh), target  , intent(inout)   :: self
      integer                 , intent(in)      :: NNew(NDIM,self % no_of_elements)
      type(FTValueDictionary) , intent(in)      :: controlVariables
      !-local-variables-----------------------------------
      integer :: eID, fID, zoneID
      logical :: saveGradients, FaceComputeQdot
      logical :: analyticalJac   ! Do we need analytical Jacobian storage?
      type(IntegerDataLinkedList_t) :: elementList
      type(IntegerDataLinkedList_t) :: facesList
      type(IntegerDataLinkedList_t) :: zoneList
      integer         , allocatable :: zoneArray(:)
      integer         , allocatable :: facesArray(:)
      integer         , allocatable :: elementArray(:)
      type(Zone_t)    , pointer :: zone
      type(Element)   , pointer :: e
      type(Face)      , pointer :: f
#if (!defined(NAVIERSTOKES))
      logical, parameter            :: computeGradients = .true.
#endif
      !---------------------------------------------------

!     **************************************
!     Check if resulting mesh is anisotropic
!     **************************************
      if ( maxval(NNew) /= minval(NNew) ) self % anisotropic = .TRUE.

      self % NDOF = 0
      do eID=1, self % no_of_elements
         self % NDOF = self % NDOF + product( NNew(:,eID) + 1 )
      end do

!     ********************
!     Some initializations
!     ********************
      saveGradients = controlVariables % logicalValueForKey("save gradients with solution")
      FaceComputeQdot = controlVariables % containsKey("acoustic analogy")

      facesList      = IntegerDataLinkedList_t(.FALSE.)
      elementList    = IntegerDataLinkedList_t(.FALSE.)
      zoneList = IntegerDataLinkedList_t(.FALSE.)
      analyticalJac  = self % storage % anJacobian

!     *********************************************
!     Adapt individual elements (geometry excluded)
!     *********************************************
!$omp parallel do schedule(runtime) private(fID, e)
      do eID=1, self % no_of_elements
         e => self % elements(eID)   ! Associate fails(!) here
         if ( all( e % Nxyz == NNew(:,eID)) ) then
            cycle
         else
            call e % pAdapt ( NNew(:,eID), self % nodeType, saveGradients, self % storage % prevSol_num )
!$omp critical
            self % Nx(eID) = NNew(1,eID)
            self % Ny(eID) = NNew(2,eID)
            self % Nz(eID) = NNew(3,eID)
            call elementList % add (eID)
            do fID=1, 6
               call facesList   % add (e % faceIDs(fID))
               if (self % faces(e % faceIDs(fID)) % FaceType  /= HMESH_BOUNDARY) then
                  call elementList % add ( mpi_partition % global2localeID (e % Connection(fID) % globID) )
               end if
            end do
!$omp end critical

         end if
!~         end associate
      end do
!$omp end parallel do

      call facesList % ExportToArray(facesArray, .TRUE.)

!     *************************
!     Adapt corresponding faces
!     *************************

!     Destruct faces storage
!     ----------------------
!$omp parallel do schedule(runtime)
      do fID=1, size(facesArray)
         call self % faces( facesArray(fID) ) % storage % destruct
      end do
!$omp end parallel do

!     Set connectivities
!     ------------------
      call self % SetConnectivitiesAndLinkFaces (self % nodeType, facesArray)

!     Construct faces storage
!     -----------------------
!$omp parallel do private(f) schedule(runtime)
      do fID=1, size(facesArray)
         f => self % faces( facesArray(fID) )  ! associate fails here in intel compilers
         call f % storage(1) % Construct(NDIM, f % Nf, f % NelLeft , computeGradients, analyticalJac, FaceComputeQdot)
         call f % storage(2) % Construct(NDIM, f % Nf, f % NelRight, computeGradients, analyticalJac, FaceComputeQdot)
      end do
!$omp end parallel do

!     ********************
!     Reconstruct geometry
!     ********************

      !* 1. Adapted elements
      !* 2. Surrounding faces of adapted elements
      !* 3. Neighbor elements of adapted elements whose intermediate face's geometry was adapted
      !* 4. Faces and elements that share a boundary with a reconstructed face (3D non-conforming representations)


      if (self % anisotropic .and. (.not. self % meshIs2D) ) then

!        Check if any of the faces belongs to a boundary
!        -----------------------------------------------
         do fID=1, size(facesArray)
            associate (f => self % faces( facesArray(fID) ) )
            if ( f % FaceType == HMESH_BOUNDARY ) then
               call zoneList % add (f % zone)
            end if
            end associate
         end do

!        Add the corresponding faces and elements
!        ----------------------------------------
         call zoneList % ExportToArray (zoneArray)

         do zoneID=1, size(zoneArray)
            zone => self % zones( zoneArray(zoneID) )    ! Compiler bug(?): If zone was implemented as associate, gfortran would not compile
            do fID=1, zone % no_of_faces
               call facesList   % add ( zone % faces(fID) )
               call elementList % add ( self % faces(zone % faces(fID)) % elementIDs(1) )
            end do
         end do
         deallocate (zoneArray   )
      end if

      deallocate ( facesArray )

      call facesList   % ExportToArray(facesArray  , .TRUE.)
      call elementList % ExportToArray(elementArray, .TRUE.)

!     Destruct old
!     ------------
      do eID=1, size (elementArray)
         call self % elements (elementArray(eID)) % geom % destruct
      end do
      do fID=1, size (facesArray)
         call self % faces (facesArray(fID)) % geom % destruct
      end do

      call self % ConstructGeometry(facesArray, elementArray)

#if defined(NAVIERSTOKES)
      call self % ComputeWallDistances(facesArray, elementArray)
#endif

!     *********
!     Finish up
!     *********
      call self % PrepareForIO

      call facesList    % destruct
      call elementList  % destruct
      call zoneList     % destruct
      nullify (zone)
      nullify (e)
      nullify (f)
      deallocate (facesArray  )
      deallocate (elementArray)

   end subroutine HexMesh_pAdapt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  HexMesh_Assign:
!  Subroutine to assign a HexMesh to another.
!  It turns out that an "impure" procedure with explicit OMP is more efficient than pure or elemental in this case.
!
   subroutine HexMesh_Assign (to, from)
      implicit none
      !-arguments----------------------------------------
      class(HexMesh), intent(inout), target :: to
      type (HexMesh), intent(in)    :: from
      !-local-variables----------------------------------
      integer :: eID
      !--------------------------------------------------
      to % numberOfFaces      = from % numberOfFaces
      to % nodeType           = from % nodeType
      to % no_of_elements     = from % no_of_elements
      to % no_of_allElements  = from % no_of_allElements
      to % dt_restriction     = from % dt_restriction
      to % NDOF               = from % NDOF

      to % meshFileName       = from % meshFileName
      to % meshIs2D           = from % meshIs2D
      to % dir2D              = from % dir2D

      to % anisotropic        = from % anisotropic
      to % child              = from % child
      to % ignoreBCnonConformities = from % child

      safedeallocate (to % Nx)
      allocate ( to % Nx ( size(from % Nx) ) )
      to % Nx = from % Nx

      safedeallocate (to % Ny)
      allocate ( to % Ny ( size(from % Ny) ) )
      to % Ny = from % Ny

      safedeallocate (to % Nz)
      allocate ( to % Nz ( size(from % Nz) ) )
      to % Nz = from % Nz

      to % storage = from % storage

      safedeallocate (to % nodes)
      allocate ( to % nodes ( size(from % nodes) ) )
!$omp parallel do schedule(runtime)
      do eID=1, size(from % nodes)
         to % nodes(eID) = from % nodes(eID)
      end do
!$omp end parallel do

      safedeallocate (to % faces)
      allocate ( to % faces ( size(from % faces) ) )
!$omp parallel do schedule(runtime)
      do eID=1, size(from % faces)
         to % faces(eID) = from % faces(eID)
      end do
!$omp end parallel do
      safedeallocate (to % elements)
      allocate ( to % elements ( size(from % elements) ) )
!$omp parallel do schedule(runtime)
      do eID=1, from % no_of_elements
         to % elements(eID) = from % elements(eID)
      end do
!$omp end parallel do

      to % MPIfaces = from % MPIfaces

      safedeallocate (to % zones)
      allocate ( to % zones ( size(from % zones) ) )
      to % zones = from % zones

!
!     Point elements' storage
!     -----------------------
!$omp parallel do schedule(runtime)
      do eID = 1, to % no_of_elements
         to % elements(eID) % storage => to % storage % elements(eID)
      end do
!$omp end parallel do

      safedeallocate(to % elements_sequential)
      safedeallocate(to % elements_mpi       )
      safedeallocate(to % faces_interior     )
      safedeallocate(to % faces_mpi          )
      safedeallocate(to % faces_boundary     )

      allocate(to % elements_sequential(size(from % elements_sequential)))
      allocate(to % elements_mpi       (size(from % elements_mpi       )))
      allocate(to % faces_interior     (size(from % faces_interior     )))
      allocate(to % faces_mpi          (size(to % faces_mpi            )))
      allocate(to % faces_boundary     (size(from % faces_boundary     )))

      to % elements_sequential = from % elements_sequential
      to % elements_mpi        = from % elements_mpi
      to % faces_interior      = from % faces_interior
      to % faces_mpi           = to % faces_mpi
      to % faces_boundary      = from % faces_boundary


   end subroutine HexMesh_Assign
END MODULE HexMeshClass
