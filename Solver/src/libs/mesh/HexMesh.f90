!
!////////////////////////////////////////////////////////////////////////
!
!      HexMesh.f95
!      Created: 2007-03-22 17:05:00 -0400 
!      By: David Kopriva  
!
!
!////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
MODULE HexMeshClass
      use SMConstants
      USE MeshTypes
      USE NodeClass
      USE ElementClass
      USE FaceClass
      use FacePatchClass
      USE TransfiniteMapClass
      use SharedBCModule
      use ElementConnectivityDefinitions
      use ZoneClass
      use PhysicsStorage
      use NodalStorageClass
      use MPI_Process_Info
      use MPI_Face_Class
      use FluidData
#if defined(NAVIERSTOKES)
      use WallDistance
#endif
#ifdef _HAS_MPI_
      use mpi
#endif
      IMPLICIT NONE

      private
      public      HexMesh, Neighbour, SurfInfo_t

      public      GetOriginalNumberOfFaces
      public      ConstructFaces, ConstructPeriodicFaces
      public      DeletePeriodicMinusFaces, GetElementsFaceIDs
!
!     ---------------
!     Mesh definition
!     ---------------
!
      type HexMesh
         integer                                   :: dir2D
         integer                                   :: numberOfFaces
         integer                                   :: nodeType
         integer                                   :: no_of_elements
         integer                                   :: no_of_allElements
         integer                                   :: dt_restriction       ! Time step restriction of last step (DT_FIXED, DT_DIFF or DT_CONV)
         character(len=LINE_LENGTH)                :: meshFileName
         type(Node)   , dimension(:), allocatable  :: nodes
         type(Face)   , dimension(:), allocatable  :: faces
         type(Element), dimension(:), allocatable  :: elements
         class(Zone_t), dimension(:), allocatable  :: zones
         logical                                   :: child       = .FALSE.         ! Is this a (multigrid) child mesh? default .FALSE.
         logical                                   :: meshIs2D    = .FALSE.         ! Is this a 2D mesh? default .FALSE.
         logical                                   :: anisotropic = .FALSE.         ! Is the mesh composed by elements with anisotropic polynomial orders? default false
         contains
            procedure :: destruct                      => DestructMesh
            procedure :: Describe                      => DescribeMesh
            procedure :: DescribePartition             => DescribeMeshPartition
            procedure :: ConstructZones                => HexMesh_ConstructZones
            procedure :: DefineAsBoundaryFaces         => HexMesh_DefineAsBoundaryFaces
            procedure :: CorrectOrderFor2DMesh         => HexMesh_CorrectOrderFor2DMesh
            procedure :: SetConnectivitiesAndLinkFaces => HexMesh_SetConnectivitiesAndLinkFaces
            procedure :: UpdateFacesWithPartition      => HexMesh_UpdateFacesWithPartition
            procedure :: ConstructGeometry             => HexMesh_ConstructGeometry
            procedure :: ProlongSolutionToFaces        => HexMesh_ProlongSolutionToFaces
            procedure :: ProlongGradientsToFaces       => HexMesh_ProlongGradientsToFaces
            procedure :: PrepareForIO                  => HexMesh_PrepareForIO
            procedure :: Export                        => HexMesh_Export
            procedure :: ExportOrders                  => HexMesh_ExportOrders
            procedure :: SaveSolution                  => HexMesh_SaveSolution
            procedure :: SaveStatistics                => HexMesh_SaveStatistics
            procedure :: ResetStatistics               => HexMesh_ResetStatistics
            procedure :: LoadSolution                  => HexMesh_LoadSolution
            procedure :: WriteCoordFile
            procedure :: UpdateMPIFacesSolution        => HexMesh_UpdateMPIFacesSolution
            procedure :: UpdateMPIFacesGradients       => HexMesh_UpdateMPIFacesGradients
            procedure :: GatherMPIFacesSolution        => HexMesh_GatherMPIFacesSolution
            procedure :: GatherMPIFacesGradients       => HexMesh_GatherMPIFacesGradients
            procedure :: FindPointWithCoords           => HexMesh_FindPointWithCoords
            procedure :: ComputeWallDistances          => HexMesh_ComputeWallDistances
      end type HexMesh

      TYPE Neighbour         ! added to introduce colored computation of numerical Jacobian (is this the best place to define this type??) - only usable for conforming meshes
         INTEGER :: elmnt(7) ! "7" hardcoded for 3D hexahedrals in conforming meshes... This definition must change if the code is expected to be more general
      END TYPE Neighbour

!
!     -------------------------------------------------------------
!     Type containing the information of the surfaces of an element
!     -> This is used only for constructung the mesh by the readers
!     -------------------------------------------------------------
      type SurfInfo_t
         ! Variables to specify that the element is a hex8
         logical         :: IsHex8 = .FALSE.
         real(kind=RP)   :: corners(NDIM,NODES_PER_ELEMENT)
         ! Variables for general elements
         type(FacePatch) :: facePatches(6)
      end type SurfInfo_t
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
      SUBROUTINE DestructMesh( self )
         IMPLICIT NONE 
         CLASS(HexMesh) :: self
         INTEGER        :: j
!
!        -----
!        Nodes
!        -----
!
         DO j = 1, SIZE( self % nodes )
            CALL DestructNode( self % nodes(j)) 
         END DO  
         DEALLOCATE( self % nodes )
!
!        --------
!        Elements
!        --------
!
         DO j = 1, SIZE(self % elements) 
            CALL DestructElement( self % elements(j) )
         END DO
         DEALLOCATE( self % elements )
!
!        -----
!        Faces
!        -----
!
         DO j = 1, SIZE(self % faces) 
            call self % faces(j) % Destruct
         END DO
         DEALLOCATE( self % faces )
!
!        -----
!        Zones
!        -----
!
         if (allocated(self % zones)) DEALLOCATE( self % zones )
         
      END SUBROUTINE DestructMesh
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
         USE FTMultiIndexTableClass 
         USE FTValueClass
         
         IMPLICIT NONE  
         TYPE(HexMesh) :: self
         LOGICAL       :: success
         
         INTEGER                 :: eID, faceNumber
         INTEGER                 :: faceID
         INTEGER                 :: nodeIDs(8), faceNodeIDs(4), j
         
         CLASS(FTMultiIndexTable), POINTER  :: table
         CLASS(FTObject), POINTER :: obj
         CLASS(FTValue) , POINTER :: v
         
         ALLOCATE(table)
         CALL table % initWithSize( SIZE( self % nodes) )
         
         self % numberOfFaces = 0
         DO eID = 1, SIZE( self % elements )
         
            nodeIDs = self % elements(eID) % nodeIDs
            DO faceNumber = 1, 6
               DO j = 1, 4
                  faceNodeIDs(j) = nodeIDs(localFaceNode(j,faceNumber)) 
               END DO
            
               IF ( table % containsKeys(faceNodeIDs) )     THEN
!
!                 --------------------------------------------------------------
!                 Add this element to the slave side of the face associated with
!                 these nodes.
!                 --------------------------------------------------------------
!
                  obj => table % objectForKeys(faceNodeIDs)
                  v   => valueFromObject(obj)
                  faceID = v % integerValue()
                  
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
          
                     CALL release(table)
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
                  ALLOCATE(v)
                  CALL v % initWithValue(self % numberOfFaces)
                  obj => v
                  CALL table % addObjectForKeys(obj,faceNodeIDs)
                  CALL release(v)
               END IF 
            END DO 
              
         END DO  
         
         CALL release(table)
         
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
      SUBROUTINE ConstructPeriodicFaces(self) 
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
      TYPE(HexMesh) :: self

! 
!-------------------- 
! Local variables 
!-------------------- 
! 
!
      REAL(KIND=RP) :: x1(NDIM), x2(NDIM)
      LOGICAL       :: master_matched(4), slave_matched(4), success
      INTEGER       :: coord, slaveNodeIDs(4), localCoord
      
      INTEGER       :: i,j,k,l 
      integer       :: zIDplus, zIDMinus, iFace, jFace
!
!     ---------------------------------------------
!     Loop to find faces with the label "periodic+"
!     ---------------------------------------------
!
!     ------------------------------
!     Loop zones with BC "periodic+"
!     ------------------------------
!
      if ( bcTypeDictionary % COUNT() .eq. 0 ) return
      do zIDPlus = 1, size(self % zones)
!
!        Cycle if the zone is not periodic+
!        ----------------------------------
         if ( trim(bcTypeDictionary % stringValueForKey(key = self % zones(zIDPlus) % Name, &
                                                      requestedLength = BC_STRING_LENGTH)) .ne. "periodic+") cycle
!
!        Reset the coordinate (changes when changing zones)
!        --------------------------------------------------
         coord = 0 
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
!           Loop zones with BC "periodic-"
!           ------------------------------
            do zIDMinus = 1, size(self % zones)
!
!              Cycle if the zone is not periodic-
!              ----------------------------------
               if ( trim(bcTypeDictionary % stringValueForKey(key = self % zones(zIDMinus) % Name, &
                                                         requestedLength = BC_STRING_LENGTH)) .ne. "periodic-") cycle
!
!              Loop faces in the periodic- zone
!              --------------------------------
mloop:         do jFace = 1, self % zones(zIDMinus) % no_of_faces
!
!                 Get the face ID
!                 ---------------
                  j = self % zones(zIDMinus) % faces(jFace)
!
!                 Consider only HMESH_UNDEFINED faces
!                 -----------------------------------
                  if ( (self % faces(j) % faceType .ne. HMESH_UNDEFINED)) cycle mloop
!
!                 ----------------------------------------------------------------------------------------
!                 The index i is a periodic+ face
!                 The index j is a periodic- face
!                 We are looking for couples of periodic+ and periodic- faces where 2 of the 3 coordinates
!                 in all the corners are shared. The non-shared coordinate has to be always the same one.
!                 ---------------------------------------------------------------------------------------
!
                  master_matched(:)   = .FALSE.     ! True if the master corner finds a partner
                  slave_matched(:)    = .FALSE.     ! True if the slave corner finds a partner
   
                  if ( coord .eq. 0 ) then
!
!                    Check all coordinates
!                    ---------------------
                     do localCoord = 1, 3
                        master_matched = .false.
                        slave_matched = .false.
mastercoord:            DO k = 1, 4
                           x1 = self%nodes(self%faces(i)%nodeIDs(k))%x                           
slavecoord:                DO l = 1, 4
                              IF (.NOT.slave_matched(l)) THEN 
                                 x2 = self%nodes(self%faces(j)%nodeIDs(l))%x        
                                 CALL CompareTwoNodes(x1, x2, master_matched(k), localCoord) 
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
!                    Check only the shared coordinates
!                    ---------------------------------
                     DO k = 1, 4
                        x1 = self%nodes(self%faces(i)%nodeIDs(k))%x                           
                        DO l = 1, 4
                           IF (.NOT.slave_matched(l)) THEN 
                              x2 = self%nodes(self%faces(j)%nodeIDs(l))%x        
                              CALL CompareTwoNodes(x1, x2, master_matched(k), coord) 
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
                           call compareTwoNodes(x1, x2, success, coord)
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
            end do            ! periodic- zones
!
!           If the code arrives here, the periodic+ face was not able to find a partner
!           ---------------------------------------------------------------------------
            print*, "When constructing periodic boundary conditions,"
            write(STD_OUT,'(A,I0,A,I0,A)') "Face ",i," in zone ",zIDPlus, &
                  " was not able to find a partner."
            errorMessage(STD_OUT)
            stop

            end do   ploop    ! periodic+ faces
         end do               ! periodic+ zones
           
      END SUBROUTINE ConstructPeriodicFaces
! 
!//////////////////////////////////////////////////////////////////////// 
!     
      SUBROUTINE CompareTwoNodes(x1, x2, success, coord) 
      use Utilities, only: almostEqual
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
      
         
      newFaceID = -1
      iFace = 0
      ALLOCATE( dummy_faces(self % numberOfFaces) )
      DO i = 1, self%numberOfFaces 
         bName = trim(bcTypeDictionary % stringValueForKey(key = self % faces(i) % boundaryName, &
                                                           requestedLength = LINE_LENGTH))
   
         if ( (self % faces(i) % faceType .ne. HMESH_UNDEFINED) .or. (trim(bName) .ne. "periodic-") ) then
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
      if ( (MPI_Process % doMPIAction) .and. MPI_Faces_Constructed ) then
         do domain = 1, MPI_Process % nProcs
            do iFace = 1, mpi_faces(domain) % no_of_faces
               mpi_faces(domain) % faceIDs(iFace) = newFaceID(mpi_faces(domain) % faceIDs(iFace))
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
      subroutine HexMesh_ProlongSolutionToFaces(self)
         implicit none
         class(HexMesh),   intent(inout)  :: self
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
            call self % elements(eID) % ProlongSolutionToFaces(self % faces(fIDs(1)),&
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
      subroutine HexMesh_ProlongGradientsToFaces(self)
         implicit none
         class(HexMesh),   intent(inout)  :: self
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
            call self % elements(eID) % ProlongGradientsToFaces(self % faces(fIDs(1)),&
                                                                self % faces(fIDs(2)),&
                                                                self % faces(fIDs(3)),&
                                                                self % faces(fIDs(4)),&
                                                                self % faces(fIDs(5)),&
                                                                self % faces(fIDs(6)) )
         end do
!$omp end do

      end subroutine HexMesh_ProlongGradientsToFaces
! 
!//////////////////////////////////////////////////////////////////////// 
! 
      subroutine HexMesh_UpdateMPIFacesSolution(self)
         use MPI_Face_Class
         implicit none
         class(HexMesh)    :: self
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
            call mpi_faces(domain) % RecvQ(domain)
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
            if ( mpi_faces(domain) % no_of_faces .eq. 0 ) cycle

            do mpifID = 1, mpi_faces(domain) % no_of_faces
               fID = mpi_faces(domain) % faceIDs(mpifID)
               thisSide = mpi_faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  mpi_faces(domain) % Qsend(counter:counter+NCONS-1) = f % storage(thisSide) % Q(:,i,j)
                  counter = counter + NCONS
               end do               ; end do
               end associate
            end do
!
!           -------------
!           Send solution
!           -------------
!
            call mpi_faces(domain) % SendQ(domain)
         end do
#endif
      end subroutine HexMesh_UpdateMPIFacesSolution

      subroutine HexMesh_UpdateMPIFacesGradients(self)
         use MPI_Face_Class
         implicit none
         class(HexMesh)    :: self
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
            call mpi_faces(domain) % RecvU_xyz(domain)
         end do
!
!        ***************
!        Gather gradients
!        ***************
!
         do domain = 1, MPI_Process % nProcs
            if ( mpi_faces(domain) % no_of_faces .eq. 0 ) cycle

            counter = 1

            do mpifID = 1, mpi_faces(domain) % no_of_faces
               fID = mpi_faces(domain) % faceIDs(mpifID)
               thisSide = mpi_faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  mpi_faces(domain) % U_xyzsend(counter:counter+NGRAD-1) = f % storage(thisSide) % U_x(:,i,j)
                  counter = counter + NGRAD
               end do               ; end do

               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  mpi_faces(domain) % U_xyzsend(counter:counter+NGRAD-1) = f % storage(thisSide) % U_y(:,i,j)
                  counter = counter + NGRAD
               end do               ; end do

               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  mpi_faces(domain) % U_xyzsend(counter:counter+NGRAD-1) = f % storage(thisSide) % U_z(:,i,j)
                  counter = counter + NGRAD
               end do               ; end do
               end associate
            end do

            call mpi_faces(domain) % SendU_xyz(domain)
         end do
#endif
      end subroutine HexMesh_UpdateMPIFacesGradients

      subroutine HexMesh_GatherMPIFacesSolution(self)
         implicit none
         class(HexMesh)    :: self
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
            call mpi_faces(domain) % WaitForSolution

            counter = 1
            do mpifID = 1, mpi_faces(domain) % no_of_faces
               fID = mpi_faces(domain) % faceIDs(mpifID)
               thisSide = mpi_faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  f % storage(otherSide(thisSide)) % Q(:,i,j) = mpi_faces(domain) % Qrecv(counter:counter+NCONS-1)
                  counter = counter + NCONS
               end do               ; end do
               end associate
            end do
         end do
#endif
      end subroutine HexMesh_GatherMPIFacesSolution

      subroutine HexMesh_GatherMPIFacesGradients(self)
         implicit none
         class(HexMesh)    :: self
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
            call mpi_faces(domain) % WaitForGradients

            counter = 1
            do mpifID = 1, mpi_faces(domain) % no_of_faces
               fID = mpi_faces(domain) % faceIDs(mpifID)
               thisSide = mpi_faces(domain) % elementSide(mpifID)
               associate(f => self % faces(fID))
               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  f % storage(otherSide(thisSide)) % U_x(:,i,j) = mpi_faces(domain) % U_xyzrecv(counter:counter+NGRAD-1)
                  counter = counter + NGRAD
               end do               ; end do

               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  f % storage(otherSide(thisSide)) % U_y(:,i,j) = mpi_faces(domain) % U_xyzrecv(counter:counter+NGRAD-1)
                  counter = counter + NGRAD
               end do               ; end do

               do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                  f % storage(otherSide(thisSide)) % U_z(:,i,j) = mpi_faces(domain) % U_xyzrecv(counter:counter+NGRAD-1)
                  counter = counter + NGRAD
               end do               ; end do
               end associate
            end do
         end do
#endif
      end subroutine HexMesh_GatherMPIFacesGradients

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
!        Get the total number of elements
!        --------------------------------
         if ( (MPI_Process % doMPIAction)) then
#ifdef _HAS_MPI_
            call mpi_allreduce(self % no_of_elements, self % no_of_allelements, &
                               1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
         else
            self % no_of_allElements = self % no_of_elements
         
         end if
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
      SUBROUTINE DescribeMesh( self , fileName, bFaceOrder )
      USE Headers
      IMPLICIT NONE
!
!--------------------------------------------------------------------
!  This subroutine describes the loaded mesh
!--------------------------------------------------------------------
!
!
!     ------------------
!     External variables
!     ------------------
!
      CLASS(HexMesh)      :: self
      CHARACTER(LEN=*)    :: fileName
      integer, intent(in) :: bFaceOrder
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER           :: fID, zoneID
      INTEGER           :: no_of_bdryfaces
      
      no_of_bdryfaces = 0
      
      write(STD_OUT,'(/)')
      call Section_Header("Reading mesh")
      write(STD_OUT,'(/)')
      
      call SubSection_Header('Mesh file "' // trim(fileName) // '".')

      write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of nodes: " , size ( self % nodes )
      write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of elements: " , size ( self % elements )
      write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of faces: " , size ( self % faces )
      
      do fID = 1 , size ( self % faces )
         if ( self % faces(fID) % faceType .ne. HMESH_INTERIOR) then
            no_of_bdryfaces = no_of_bdryfaces + 1
         end if
      end do

      write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of boundary faces: " , no_of_bdryfaces
      write(STD_OUT,'(30X,A,A28,I10)') "->" , "Order of curved faces: " , bFaceOrder
      
!     Describe the zones
!     ------------------
      write(STD_OUT,'(/)')
      call Section_Header("Creating zones")
      write(STD_OUT,'(/)')
      
      do zoneID = 1, size(self % zones)
         WRITE(STD_OUT,'(30X,A,A7,I0,A15,A)') "->", "  Zone ",zoneID, " for boundary: ",trim(self % zones(zoneID) % Name)
         write(STD_OUT,'(32X,A28,I0)') 'Number of faces: ', self % zones(zoneID) % no_of_faces
      end do
      
      END SUBROUTINE DescribeMesh     

      SUBROUTINE DescribeMeshPartition( self , fileName )
      USE Headers
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
      CHARACTER(LEN=*)  :: fileName
#ifdef _HAS_MPI_
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER           :: fID, zoneID, rank, ierr
      INTEGER           :: no_of_bdryfaces, no_of_mpifaces
      integer           :: no_of_nodesP(MPI_Process % nProcs)
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
      call mpi_gather(size(self % nodes)    , 1 , MPI_INT , no_of_nodesP    , 1 , MPI_INT , 0 , MPI_COMM_WORLD , ierr)
      call mpi_gather(size(self % faces)    , 1 , MPI_INT , no_of_facesP    , 1 , MPI_INT , 0 , MPI_COMM_WORLD , ierr)
      call mpi_gather(no_of_bdryfaces       , 1 , MPI_INT , no_of_bfacesP   , 1 , MPI_INT , 0 , MPI_COMM_WORLD , ierr)
      call mpi_gather(no_of_mpifaces        , 1 , MPI_INT , no_of_mpifacesP , 1 , MPI_INT , 0 , MPI_COMM_WORLD , ierr)

      if ( .not. MPI_Process % isRoot ) return
       
      write(STD_OUT,'(/)')
      call Section_Header("Mesh partitions")
      write(STD_OUT,'(/)')
      
      do rank = 1, MPI_Process % nProcs 

         write(partitionID,'(A,I0)') "Partition ", rank
         call SubSection_Header(trim(partitionID))

         write(STD_OUT,'(30X,A,A28,I10)') "->" , "Number of nodes: " , no_of_nodesP(rank)
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
      SUBROUTINE WriteCoordFile(self,FileName)
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
         CHARACTER(len=*)     :: FileName    !<  ...
         !--------------------------------------------------------
         INTEGER              :: NumOfElem
         INTEGER              :: i, j, k, el, Nx, Ny, Nz, ndof, cooh
         !--------------------------------------------------------
          
         NumOfElem = SIZE(self % elements)
!
!        ------------------------------------------------------------------------
!        Determine the number of degrees of freedom
!           TODO: Move this to another place if needed in other parts of the code
!        ------------------------------------------------------------------------
!
         ndof = 0
         DO el = 1, NumOfElem
            Nx = self % elements(el) % Nxyz(1)
            Ny = self % elements(el) % Nxyz(2)
            Nz = self % elements(el) % Nxyz(3)
            ndof = ndof + (Nx + 1)*(Ny + 1)*(Nz + 1)*NCONS
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
!        This subroutine gets the 2D direction in the local frame for each
!     element and sets the polynomial order to zero in that direction
!     --------------------------------------------------------------------
!
!//////////////////////////////////////////////////////////////////////////////
!
      subroutine HexMesh_CorrectOrderFor2DMesh(self, dir2D)
         use Utilities, only: almostEqual
         implicit none
         class(HexMesh),   intent(inout) :: self
         integer,          intent(in)    :: dir2D
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: eID, nID, no_of_orientedNodes
         integer  :: face1Nodes(NODES_PER_FACE)
         integer  :: face2Nodes(NODES_PER_FACE)
         real(kind=RP)  :: d2D(NDIM)
         real(kind=RP)  :: xNodesF1(NDIM,NODES_PER_FACE)
         real(kind=RP)  :: xNodesF2(NDIM,NODES_PER_FACE)
         real(kind=RP)  :: dx(NDIM,NODES_PER_FACE)

         self % dir2D = dir2D
!
!        Construct the normal vector
!        ---------------------------
         select case (dir2D)
         case(IX)
            d2D = [1.0_RP, 0.0_RP, 0.0_RP]

         case(IY)
            d2D = [0.0_RP, 1.0_RP, 0.0_RP]

         case(IZ)
            d2D = [0.0_RP, 0.0_RP, 1.0_RP]

         end select

         do eID = 1, self % no_of_elements
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
!           Check how many delta x vectors are parallel to the 2D direction
!           ---------------------------------------------------------------
            no_of_orientedNodes = 0
            do nID = 1, NODES_PER_FACE
               if ( almostEqual(abs(dot_product(dx(:,nID),d2D)),1.0_RP) ) then
                  no_of_orientedNodes = no_of_orientedNodes + 1 
               elseif ( almostEqual(abs(dot_product(dx(:,nID),d2D)),0.0_RP)) then
                  no_of_orientedNodes = no_of_orientedNodes - 1
               end if
            end do
!
!           Check if the direction is the geometrical 2D direction
!           ------------------------------------------------------
            if ( no_of_orientedNodes .eq. 4 ) then
!
!              This is the 2D direction
!              ------------------------
               e % dir2D = IX
               e % Nxyz(1) = 0
               e % spAxi => NodalStorage(0)
            
            elseif ( no_of_orientedNodes .ne. -4 ) then
!
!              If all delta x vectors are not perpendicular to the 2D direction, the mesh is not 2D
!              ------------------------------------------------------------------------------------
               print*, "The mesh does not seem to be 2D for the given direction"
               errorMessage(STD_OUT)
               stop

            end if
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
!           Check how many delta x vectors are parallel to the 2D direction
!           ---------------------------------------------------------------
            no_of_orientedNodes = 0
            do nID = 1, NODES_PER_FACE
               if ( almostEqual(abs(dot_product(dx(:,nID),d2D)),1.0_RP) ) then
                  no_of_orientedNodes = no_of_orientedNodes + 1 
               elseif ( almostEqual(abs(dot_product(dx(:,nID),d2D)),0.0_RP)) then
                  no_of_orientedNodes = no_of_orientedNodes - 1
               end if
            end do
!
!           Check if the direction is the geometrical 2D direction
!           ------------------------------------------------------
            if ( no_of_orientedNodes .eq. 4 ) then
!
!              This is the 2D direction
!              ------------------------
               e % dir2D = IY
               e % Nxyz(2) = 0
               e % spAeta => NodalStorage(0)
            
            elseif ( no_of_orientedNodes .ne. -4 ) then
!
!              If all delta x vectors are not perpendicular to the 2D direction, the mesh is not 2D
!              ------------------------------------------------------------------------------------
               print*, "The mesh does not seem to be 2D for the given direction"
               errorMessage(STD_OUT)
               stop

            end if
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
!           Check how many delta x vectors are parallel to the 2D direction
!           ---------------------------------------------------------------
            no_of_orientedNodes = 0
            do nID = 1, NODES_PER_FACE
               if ( almostEqual(abs(dot_product(dx(:,nID),d2D)),1.0_RP) ) then
                  no_of_orientedNodes = no_of_orientedNodes + 1 
               elseif ( almostEqual(abs(dot_product(dx(:,nID),d2D)),0.0_RP)) then
                  no_of_orientedNodes = no_of_orientedNodes - 1
               end if
            end do
!
!           Check if the direction is the geometrical 2D direction
!           ------------------------------------------------------
            if ( no_of_orientedNodes .eq. 4 ) then
!
!              This is the 2D direction
!              ------------------------
               e % dir2D = IZ
               e % Nxyz(3) = 0
               e % spAzeta => NodalStorage(0)
            
            elseif ( no_of_orientedNodes .ne. -4 ) then
!
!              If all delta x vectors are not perpendicular to the 2D direction, the mesh is not 2D
!              ------------------------------------------------------------------------------------
               print*, "The mesh does not seem to be 2D for the given direction"
               errorMessage(STD_OUT)
               stop

            end if
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
      subroutine HexMesh_SetConnectivitiesAndLinkFaces(self,nodes)
         implicit none
         class(HexMesh)       :: self
         integer              :: nodes
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: fID, SideL, SideR
         integer  :: NelL(2), NelR(2), k ! Polynomial orders on left and right of a face
         integer  :: domain, MPI_NDOFS(MPI_Process % nProcs), mpifID

         do fID = 1, size(self % faces)
            associate(f => self % faces(fID))
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
               call eL % Connection(SideL) % Construct(1)
               eL % Connection( SideL ) % ElementIDs(1) = eR % eID

               eR % NumberOfConnections(SideR) = 1
               call eR % Connection(SideR) % Construct(1)
               eR % Connection( SideR ) % ElementIDs(1) = eL % eID
               end associate

            case (HMESH_BOUNDARY)
               associate(eL => self % elements(f % elementIDs(1)))
!
!              Get polynomial orders of elements
!              ---------------------------------
               NelL = eL % Nxyz(axisMap(:, f % elementSide(1)))
               NelR = NelL
               end associate

            case (HMESH_MPI)
!
!              ****************************************************
!              TODO account for different orders accross both sides
!              ****************************************************
!
               associate(e => self % elements(maxval(f % elementIDs)))
               NelL = e % Nxyz(axisMap(:,maxval(f % elementSide)))
               NelR = NelL
               end associate

            end select
         
            call f % LinkWithElements(NCONS, NGRAD, NelL, NelR, nodes)
            
            end associate
         end do
!
!        --------------------------
!        Allocate MPI Faces storage
!        --------------------------
!
         if ( MPI_Process % doMPIAction ) then
            if ( .not. allocated(mpi_faces) ) return
#if _HAS_MPI_
            MPI_NDOFS = 0

            do domain = 1, MPI_Process % nProcs
               do mpifID = 1, mpi_faces(domain) % no_of_faces
                  fID = mpi_faces(domain) % faceIDs(mpifID)
                  associate( f => self % faces(fID) ) 
                  MPI_NDOFS(domain) = MPI_NDOFS(domain) + product(f % Nf + 1)
                  end associate
               end do
            end do

            call ConstructMPIFacesStorage(NCONS, NGRAD, MPI_NDOFS)
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
               call mpi_faces(domain) % Construct(no_of_mpifaces(domain))
            end if
         end do

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
            f % elementIDs(otherSide(eSide)) = HMESH_NONE
            f % elementSide(eSide) = side
            f % elementSide(otherSide(eSide)) = HMESH_NONE
            end associate
!
!           Create MPI Face
!           ---------------
            domain = partition % mpiface_sharedDomain(bFace)
            no_of_mpifaces(domain) = no_of_mpifaces(domain) + 1
            mpi_faces(domain) % faceIDs(no_of_mpifaces(domain)) = fID
            mpi_faces(domain) % elementSide(no_of_mpifaces(domain)) = eSide

         end do

      end subroutine HexMesh_UpdateFacesWithPartition
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------------------
!     Construct geometry of faces and elements
!     -> This routine guarantees that the mapping is subparametric or isoparametric
!     -----------------------------------------------------------------------------
      subroutine HexMesh_ConstructGeometry(self,SurfInfo)
         implicit none
         !--------------------------------
         class(HexMesh)    , intent(inout) :: self
         type(SurfInfo_t)  , intent(inout) :: SurfInfo(:)
         !--------------------------------
         integer                    :: i
         integer                    :: fID, eID                        ! Face and element counters
         integer                    :: eIDLeft, eIDRight, e            ! Element IDs on left and right of a face
         integer                    :: SideIDL, SideIDR, side          ! Side of elements on left and right of a face
         integer                    :: buffer                          ! A temporal variable
         integer                    :: NSurfL(2), NSurfR(2), NSurf(2)  ! Polynomial order the face was constructef with
         integer                    :: Nelf(2), Nel(2), rot
         integer                    :: CLN(2)                          ! Chebyshev-Lobatto face orders
         REAL(KIND=RP)              :: corners(NDIM,NODES_PER_ELEMENT) ! Variable just for initializing purposes
         real(kind=RP), allocatable :: faceCL(:,:,:)                   ! Coordinates of the Chebyshev-Lobatto nodes on the face
         
         type(TransfiniteHexMap), pointer :: hexMap, hex8Map, genHexMap
         !--------------------------------
         
         integer  :: zoneID, zonefID
         integer, allocatable :: bfOrder(:)
         corners = 0._RP

!
!        ******************************************************************
!        Find the polynomial order of the boundaries for anisotropic meshes
!        -> Unlike isotropic meshes, anisotropic meshes need boundary orders
!           bfOrder=N-1 in 3D (not in 2D) to be free-stream-preserving
!        ******************************************************************
!
         if (self % anisotropic .and. (.not. self % meshIs2D) ) then
            allocate ( bfOrder(size(self % zones)) )
            bfOrder = 50 ! Initialize to a big number
            do zoneID=1, size(self % zones)
               do zonefID = 1, self % zones(zoneID) % no_of_faces
                  fID = self % zones(zoneID) % faces(zonefID)
                  
                  associate( f => self % faces(fID) )
                  bfOrder(zoneID) = min(bfOrder(zoneID),f % NfLeft(1),f % NfLeft(2))
                  end associate
               end do
               bfOrder(zoneID) = bfOrder(zoneID) - 1
               call NodalStorage (bfOrder(zoneID)) % construct (self % nodeType, bfOrder(zoneID))
            end do
         end if
         
!
!        **************************************************************
!        Check surfaces' integrity and adapt them to the solution order
!        **************************************************************
!
         do fID=1, size(self % faces)
            associate( f => self % faces(fID) )
            
!
!           Check if the surfaces description in mesh file is consistent 
!           ------------------------------------------------------------
            select case (f % faceType)

            case (HMESH_INTERIOR)
               eIDLeft  = f % elementIDs(1)
               SideIDL  = f % elementSide(1)
               NSurfL   = SurfInfo(eIDLeft)  % facePatches(SideIDL) % noOfKnots - 1
            
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
                  write(STD_OUT,*) '   N Left:  ', SurfInfo(eIDLeft) % facePatches(SideIDL) % noOfKnots - 1
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
               case ( 1, 3, 4, 6 ) ! Local x and y axis are perpendicular  ! TODO this is correct?
                  if (CLN(1) /= CLN(2)) then
                     buffer = CLN(1)
                     CLN(1) = CLN(2)
                     CLN(2) = buffer
                  end if
               end select

               if ( any(CLN < NSurfR) ) then       ! TODO JMT: I have added this.. is correct?
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

               if     (SurfInfo(eIDLeft)  % IsHex8 .or. all(NSurfL == 1)) cycle
               
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
               side = maxval(f % elementSide)
               NSurf = SurfInfo(eID) % facePatches(side) % noOfKnots - 1
            
               if ( SurfInfo(eID) % IsHex8 .or. all(NSurf == 1) ) cycle

               CLN(1) = f % NfLeft(1)  ! TODO in MPI faces, p-adaption has
               CLN(2) = f % NfLeft(2)  ! not been accounted yet.

               if ( side .eq. 2 ) then    ! Right faces need to be rotated
                  select case ( f % rotation )
                  case ( 1, 3, 4, 6 ) ! Local x and y axis are perpendicular  ! TODO this is correct?
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
!
!        ----------------------------
!        Construct elements' geometry
!        ----------------------------
!
         allocate(hex8Map)
         call hex8Map % constructWithCorners(corners)
         allocate(genHexMap)
         
         do eID=1, size(self % elements)
            
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
         do fID=1, size(self % faces)
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
         do fID = 1, size(self % faces)
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

         if ( .not. MPI_Faces_Constructed ) return
!
!        Get the maximum number of faces
!        -------------------------------
         no_of_max_faces = 0
         do i = 1, MPI_Process % nProcs
            no_of_max_faces = max(no_of_max_faces, mpi_faces(i) % no_of_faces)
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
            if ( mpi_faces(i) % no_of_faces .le. 0 ) then
               recvReq(i) = MPI_REQUEST_NULL
               cycle

            end if

            call mpi_irecv(hrecv(:,i), mpi_faces(i) % no_of_faces, MPI_DOUBLE, i-1, MPI_ANY_TAG, &
                           MPI_COMM_WORLD, recvReq(i), ierr)

         end do

!
!        Gather the distances to send
!        ------------------------------
         do i = 1, MPI_Process % nProcs
            if ( mpi_faces(i) % no_of_faces .le. 0 ) cycle

            do mpifID = 1, mpi_faces(i) % no_of_faces
               fID = mpi_faces(i) % faceIDs(mpifID)
               hsend(mpifID,i) = self % faces(fID) % geom % h
            end do
         end do
!
!        Send the distances
!        ------------------
         do i = 1, MPI_Process % nProcs
            if ( mpi_faces(i) % no_of_faces .le. 0 ) then
               sendReq(i) = MPI_REQUEST_NULL
               cycle

            end if

            call mpi_isend(hsend(:,i), mpi_faces(i) % no_of_faces, MPI_DOUBLE, i-1, DEFAULT_TAG, &
                           MPI_COMM_WORLD, sendReq(i), ierr)

         end do

         call mpi_waitall(MPI_Process % nProcs, recvReq, recvSt, ierr)
!
!        Collect the distances
!        ---------------------      
         do i = 1, MPI_Process % nProcs
            if ( mpi_faces(i) % no_of_faces .le. 0 ) cycle

            do mpifID = 1, mpi_faces(i) % no_of_faces
               fID = mpi_faces(i) % faceIDs(mpifID)
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
         integer        :: fID, eID, pos
         character(len=LINE_LENGTH)    :: meshName
         real(kind=RP), parameter      :: refs(NO_OF_SAVED_REFS) = 0.0_RP
         interface
            character(len=LINE_LENGTH) function RemovePath( inputLine )
               use SMConstants
               implicit none
               character(len=*)     :: inputLine
            end function RemovePath
      
            character(len=LINE_LENGTH) function getFileName( inputLine )
               use SMConstants
               implicit none
               character(len=*)     :: inputLine
            end function getFileName
         end interface

            
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
            pos = POS_INIT_DATA + (e % globID-1)*5*SIZEOF_INT + 3*e % offsetIO*SIZEOF_RP
            call writeArray(fID, e % geom % x, position=pos)
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
         interface
            character(len=LINE_LENGTH) function RemovePath( inputLine )
               use SMConstants
               implicit none
               character(len=*)     :: inputLine
            end function RemovePath
      
            character(len=LINE_LENGTH) function getFileName( inputLine )
               use SMConstants
               implicit none
               character(len=*)     :: inputLine
            end function getFileName
         end interface
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
#if defined(NAVIERSTOKES)
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
         integer  :: fid, eID, pos, padding
         real(kind=RP)                    :: refs(NO_OF_SAVED_REFS) 
!
!        Gather reference quantities
!        ---------------------------
         refs(GAMMA_REF) = thermodynamics % gamma
         refs(RGAS_REF)  = thermodynamics % R
         refs(RHO_REF)   = refValues      % rho
         refs(V_REF)     = refValues      % V
         refs(T_REF)     = refValues      % T
         refs(MACH_REF)  = dimensionless  % Mach
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
            pos = POS_INIT_DATA + (e % globID-1)*5*SIZEOF_INT + padding*e % offsetIO * SIZEOF_RP
            call writeArray(fid, e % storage % Q, position=pos)
            if ( saveGradients .and. computeGradients ) then
               write(fid) e % storage % U_x
               write(fid) e % storage % U_y
               write(fid) e % storage % U_z
            end if
            end associate
         end do
         close(fid)
!
!        Close the file
!        --------------
         call SealSolutionFile(trim(name))

      end subroutine HexMesh_SaveSolution
#elif defined(CAHNHILLIARD)
!
!     **************************************************************************
!           Save solution subroutine for the Cahn-Hilliard equations. It saves
!        the concentration and the chemical potential
!     **************************************************************************
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
         integer  :: fid, eID, pos, padding
         real(kind=RP)                    :: refs(NO_OF_SAVED_REFS) 
!
!        Dummy references
!        ----------------
         refs = 0.0_RP
!
!        Create new file
!        ---------------
         call CreateNewSolutionFile(trim(name),SOLUTION_CAHNHILLIARD_FILE, self % nodeType, &
                                    self % no_of_allElements, iter, time, refs)
         padding = 1
!
!        Write arrays
!        ------------
         fID = putSolutionFileInWriteDataMode(trim(name))
         do eID = 1, self % no_of_elements
            associate( e => self % elements(eID) )
            pos = POS_INIT_DATA + (e % globID-1)*5*SIZEOF_INT + padding*e % offsetIO * SIZEOF_RP
            call writeArray(fid, e % storage % Q, position=pos)
            !write(fid) e % storage % mu
            end associate
         end do
         close(fid)
!
!        Close the file
!        --------------
         call SealSolutionFile(trim(name))

      end subroutine HexMesh_SaveSolution
#endif

      subroutine HexMesh_SaveStatistics(self, iter, time, name)
         use SolutionFile
         implicit none
         class(HexMesh),      intent(in)        :: self
         integer,             intent(in)        :: iter
         real(kind=RP),       intent(in)        :: time
         character(len=*),    intent(in)        :: name
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: fid, eID, pos
         real(kind=RP)                    :: refs(NO_OF_SAVED_REFS) 
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
#else
         refs = 0.0_RP
#endif
!
!        Create new file
!        ---------------
         call CreateNewSolutionFile(trim(name),STATS_FILE, self % nodeType, self % no_of_allElements, iter, time, refs)
!
!        Write arrays
!        ------------
         fID = putSolutionFileInWriteDataMode(trim(name))
         do eID = 1, self % no_of_elements
            associate( e => self % elements(eID) )
            pos = POS_INIT_DATA + (e % globID-1)*5*SIZEOF_INT + 9*e % offsetIO*SIZEOF_RP
            call writeArray(fid, e % storage % stats % data, position=pos)
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

      subroutine HexMesh_LoadSolution( self, fileName, initial_iteration, initial_time ) 
         use SolutionFile
         IMPLICIT NONE
         CLASS(HexMesh)             :: self
         character(len=*)           :: fileName
         integer,       intent(out) :: initial_iteration
         real(kind=RP), intent(out) :: initial_time
         
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER          :: fID, eID, fileType, no_of_elements, flag, nodetype
         integer          :: padding, pos
         integer          :: Nxp1, Nyp1, Nzp1, no_of_eqs, array_rank
         character(len=SOLFILE_STR_LEN)      :: rstName
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
            padding = NCONS + 3 * NGRAD

         case(STATS_FILE)
            print*, "The selected restart file is a statistics file"
            errorMessage(STD_OUT)
            stop
         case(SOLUTION_CAHNHILLIARD_FILE)
            padding = 1*NCONS
         case default
            print*, "Unknown restart file format"
            errorMessage(STD_OUT)
            stop
         end select
!
!        Get the node type
!        -----------------
         nodeType = getSolutionFileNodeType(trim(fileName))

         if ( nodeType .ne. self % nodeType ) then
            print*, "Solution file uses a different discretization nodes that the mesh."
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

            read(fID) e % storage % Q 
            end associate
         end do
!
!        Close the file
!        --------------
         close(fID)

      END SUBROUTINE HexMesh_LoadSolution
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
         logical     :: success

         HexMesh_FindPointWithCoords = .false.

         if ( present(optionalElements) ) then
            do op_eID = 1, size(optionalElements)
               if ( optionalElements(op_eID) .eq. -1 ) cycle
               associate(e => self % elements(optionalElements(op_eID)))
               success = e % FindPointWithCoords(x, xi) 
               if ( success ) then
                  eID = optionalElements(op_eID)
                  HexMesh_FindPointWithCoords = .true.
                  return 
               end if
               end associate
            end do
         end if      
         
         do eID = 1, self % no_of_elements
            associate(e => self % elements(eID))
            success = e % FindPointWithCoords(x, xi)
            if ( success ) then
               HexMesh_FindPointWithCoords = .true.
               return
            end if
            end associate
         end do

      end function HexMesh_FindPointWithCoords

      subroutine HexMesh_ComputeWallDistances(self)
         implicit none
         class(HexMesh)     :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: eID, i, j, k, no_of_wallDOFS
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
         do eID = 1, self % no_of_elements
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
         do eID = 1, size(self % faces)
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
         character(len=LINE_LENGTH) :: zoneName

         no_of_localWallDOFS = 0
         do zID = 1, size(self % zones)
            zoneName = trim(bcTypeDictionary % stringValueForKey(trim(self % zones(zID) % Name), LINE_LENGTH))

            if ( (trim(zoneName) .ne. "noslipadiabaticwall") .and. (trim(zoneName) .ne. "noslipisothermalwall") ) then
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
            zoneName = trim(bcTypeDictionary % stringValueForKey(trim(self % zones(zID) % Name), LINE_LENGTH))

            if ( (trim(zoneName) .ne. "noslipadiabaticwall") .and. (trim(zoneName) .ne. "noslipisothermalwall") ) then
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
END MODULE HexMeshClass
      
