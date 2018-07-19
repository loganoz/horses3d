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
MODULE Read_SpecMesh
      use SMConstants
      use MeshTypes
      use ElementConnectivityDefinitions
      USE TransfiniteMapClass
      use FacePatchClass
      use MappedGeometryClass
      use NodeClass
      use ElementClass
      use HexMeshClass
      use sharedBCModule
      use PhysicsStorage
      use FileReadingUtilities      , only: getFileName
      implicit none
      
      private
      public ConstructMesh_FromSpecMeshFile_, NumOfElems_SpecMesh
      
!
!     ========
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
!!    Constructs mesh from mesh file
!!    Only valid for conforming meshes
      SUBROUTINE ConstructMesh_FromSpecMeshFile_( self, fileName, nodes, Nx, Ny, Nz, dir2D, success, export )
         USE Physics
         use PartitionedMeshClass
         use MPI_Process_Info
         use Utilities, only: UnusedUnit
         IMPLICIT NONE
!
!        ---------------
!        Input variables
!        ---------------
!
         type(HexMesh)     , intent(out)  :: self
         integer           , intent(in)   :: nodes
         CHARACTER(LEN=*)  , intent(in)   :: fileName
         integer           , intent(in)   :: Nx(:), Ny(:), Nz(:)     !<  Polynomial orders for all the elements
         integer           , intent(in)   :: dir2D
         LOGICAL           , intent(out)  :: success
         logical, optional , intent(in)   :: export 
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                         :: numberOfElements
         integer                         :: numberOfNodes
         integer                         :: numberOfBoundaryFaces
         integer                         :: numberOfFaces
         
         integer                         :: bFaceOrder, numBFacePoints
         integer                         :: i, j, k, l
         integer                         :: fUnit, fileStat
         integer                         :: nodeIDs(NODES_PER_ELEMENT), nodeMap(NODES_PER_FACE)
         real(kind=RP)                   :: x(NDIM)
         integer                         :: faceFlags(FACES_PER_ELEMENT)
         CHARACTER(LEN=BC_STRING_LENGTH) :: names(FACES_PER_ELEMENT)
         TYPE(FacePatch), DIMENSION(6)   :: facePatches
         real(kind=RP)                   :: corners(NDIM,NODES_PER_ELEMENT)
         logical                         :: export_mesh
!
!        ------------------
!        For curved patches
!        ------------------
!
         type(SurfInfo_t), allocatable                  :: SurfInfo(:)
         real(kind=RP)  , DIMENSION(:)    , ALLOCATABLE :: uNodes, vNodes
         real(kind=RP)  , DIMENSION(:,:,:), ALLOCATABLE :: values
!
!        ----------------
!        For flat patches
!        ----------------
!
         real(kind=RP)  , DIMENSION(2)     :: uNodesFlat = [-1.0_RP,1.0_RP]
         real(kind=RP)  , DIMENSION(2)     :: vNodesFlat = [-1.0_RP,1.0_RP]
         real(kind=RP)  , DIMENSION(3,2,2) :: valuesFlat
!
!        ********************************
!        Check if a mesh partition exists
!        ********************************
!
         if ( mpi_partition % Constructed ) then
            call ConstructMeshPartition_FromSpecMeshFile_( self, fileName, nodes, Nx, Ny, Nz, dir2D, success ) 
            return
         end if
          
         numberOfBoundaryFaces = 0
         success               = .TRUE.
!
!        -----------------------
!        Read header information
!        -----------------------
!
         fUnit = UnusedUnit()
         OPEN( UNIT = fUnit, FILE = fileName, iostat = fileStat )
         IF ( fileStat /= 0 )     THEN
            PRINT *, "Error opening file: ", fileName
            success = .FALSE.
            RETURN 
         END IF
         
         READ(fUnit,*) numberOfNodes, numberOfElements, bFaceOrder

         self % nodeType = nodes
         self % no_of_elements = numberOfElements
!
!        ----------------------------
!        Set up for face patches
!        Face patches are defined 
!        at chebyshev- lobatto points
!        ----------------------------
!
         numBFacePoints = bFaceOrder + 1
         allocate(uNodes(numBFacePoints))
         allocate(vNodes(numBFacePoints))
         allocate(values(3,numBFacePoints,numBFacePoints))
         
         DO i = 1, numBFacePoints
            uNodes(i) = -COS((i-1.0_RP)*PI/(numBFacePoints-1.0_RP)) 
            vNodes(i) = -COS((i-1.0_RP)*PI/(numBFacePoints-1.0_RP)) 
         END DO
         
         DO k = 1, 6 ! Most patches will be flat, so set up for self
            CALL facePatches(k) % construct(uNodesFlat,vNodesFlat) 
         END DO  
!
!        ---------------
!        Allocate memory
!        ---------------
!
         allocate( self % elements(numberOfelements) )
         allocate( self % nodes(numberOfNodes) )
         allocate( SurfInfo(numberOfelements) )

!
!        ----------------------------------
!        Read nodes: Nodes have the format:
!        x y z
!        ----------------------------------
!
         DO j = 1, numberOfNodes 
            READ( fUnit, * ) x
            CALL ConstructNode( self % nodes(j), x, j )
         END DO
!
!        -----------------------------------------
!        Read elements: Elements have the format:
!        node1ID node2ID node3ID node4ID ... node8ID
!        b1 b2 b3 b4 b5 b6
!           (=0 for straight side, 1 for curved)
!        If curved boundaries, then for each:
!           for j = 0 to bFaceOrder
!              x_j  y_j  z_j
!           next j
!        bname1 bname2 bname3 bname4 bname5 bname6
!        -----------------------------------------
!
         DO l = 1, numberOfElements 
            READ( fUnit, * ) nodeIDs
            READ( fUnit, * ) faceFlags
!
!           -----------------------------------------------------------------------------
!           If the faceFlags are all zero, then self is a straight-sided
!           hex. In that case, set the corners of the hex8Map and use that in determining
!           the element geometry.
!           -----------------------------------------------------------------------------
!
            IF(MAXVAL(faceFlags) == 0)     THEN
               DO k = 1, NODES_PER_ELEMENT
                  corners(:,k) = self % nodes(nodeIDs(k)) % x
               END DO
               SurfInfo(l) % IsHex8 = .TRUE.
               SurfInfo(l) % corners = corners
               
            ELSE
!
!              --------------------------------------------------------------
!              Otherwise, we have to look at each of the faces of the element 
!              --------------------------------------------------------------
!
               DO k = 1, FACES_PER_ELEMENT 
                 IF ( faceFlags(k) == 0 )     THEN
!
!                    ----------
!                    Flat faces
!                    ----------
!
                     nodeMap           = localFaceNode(:,k)
                     valuesFlat(:,1,1) = self % nodes(nodeIDs(nodeMap(1))) % x
                     valuesFlat(:,2,1) = self % nodes(nodeIDs(nodeMap(2))) % x
                     valuesFlat(:,2,2) = self % nodes(nodeIDs(nodeMap(3))) % x
                     valuesFlat(:,1,2) = self % nodes(nodeIDs(nodeMap(4))) % x
                     
                     call SurfInfo(l) % facePatches(k) % construct(uNodesFlat, vNodesFlat, valuesFlat)
                     
                  ELSE
!
!                    -------------
!                    Curved faces 
!                    -------------
!
                     DO j = 1, numBFacePoints
                        DO i = 1, numBFacePoints
                           READ( fUnit, * ) values(:,i,j)
                        END DO  
                     END DO
                     
                     call SurfInfo(l) % facePatches(k) % construct(uNodes, vNodes, values)
                     
                  END IF
               END DO
               
            END IF 
!
!           -------------------------
!           Now construct the element
!           -------------------------
!
            call self % elements(l) % Construct (Nx(l), Ny(l), Nz(l), nodeIDs , l, l)
            
            READ( fUnit, * ) names
            CALL SetElementBoundaryNames( self % elements(l), names )
            
            DO k = 1, 6
               IF(TRIM(names(k)) /= emptyBCName) then
                  numberOfBoundaryFaces = numberOfBoundaryFaces + 1
                  if ( all(trim(names(k)) .ne. zoneNameDictionary % allKeys()) ) then
                     call zoneNameDictionary % addValueForKey(trim(names(k)), trim(names(k)))
                  end if
               end if
            END DO  
         END DO      ! l = 1, numberOfElements
        
!
!        ---------------------------
!        Construct the element faces
!        ---------------------------
!
         numberOfFaces        = (6*numberOfElements + numberOfBoundaryFaces)/2
         self % numberOfFaces = numberOfFaces
         allocate( self % faces(self % numberOfFaces) )
         CALL ConstructFaces( self, success )
!
!        -------------------------
!        Build the different zones
!        -------------------------
!
         call self % ConstructZones()
!
!        ---------------------------
!        Construct periodic faces
!        ---------------------------
!
         CALL ConstructPeriodicFaces( self )
!
!        ---------------------------
!        Delete periodic- faces
!        ---------------------------
!
         CALL DeletePeriodicMinusFaces( self )
!
!        ---------------------------
!        Assign faces ID to elements
!        ---------------------------
!
         CALL getElementsFaceIDs(self)
!
!        ---------------------
!        Define boundary faces
!        ---------------------
!
         call self % DefineAsBoundaryFaces()
!
!        -----------------------------------
!        Check if this is a 2D extruded mesh
!        -----------------------------------
!
         call self % CheckIfMeshIs2D()
!
!        -------------------------------
!        Set the mesh as 2D if requested
!        -------------------------------
!
         if ( dir2D .ne. 0 ) then
            call SetMappingsToCrossProduct
            call self % CorrectOrderFor2DMesh(dir2D)
         end if
!
!        ------------------------------
!        Set the element connectivities
!        ------------------------------
!
         call self % SetConnectivitiesAndLinkFaces(nodes)
!
!        ---------------------------------------
!        Construct elements' and faces' geometry
!        ---------------------------------------
!
         call self % ConstructGeometry(SurfInfo)
            
         CLOSE( fUnit )
!
!        ---------
!        Finish up
!        ---------
!
         if (.not. self % child) CALL self % Describe( trim(fileName), bFaceOrder )
!
!        -------------------------------------------------------------
!        Prepare mesh for I/O only if the code is running sequentially
!        -------------------------------------------------------------
!
         if (present(export)) then
            export_mesh = export
         else
            export_mesh = .true.
         end if
         if ( .not. MPI_Process % doMPIAction ) then
            call self % PrepareForIO
            if ((.not. self % child) .and. export_mesh) then
               call self % Export( trim(fileName) )
            end if
         end if
         
      END SUBROUTINE ConstructMesh_FromSpecMeshFile_

      SUBROUTINE ConstructMeshPartition_FromSpecMeshFile_( self, fileName, nodes, Nx, Ny, Nz, dir2D, success )
         USE Physics
         use PartitionedMeshClass
         use MPI_Process_Info
         use MPI_Face_Class
         use Utilities, only: UnusedUnit
         IMPLICIT NONE
!
!        ---------------
!        Input variables
!        ---------------
!
         CLASS(HexMesh)      :: self
         integer             :: nodes
         CHARACTER(LEN=*)    :: fileName
         integer             :: Nx(:), Ny(:), Nz(:)     !<  Polynomial orders for all the elements
         integer, intent(in) :: dir2D
         LOGICAL             :: success
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                         :: numberOfAllElements
         integer                         :: numberOfAllNodes
         integer                         :: numberOfFaces
         integer, allocatable            :: globalToLocalNodeID(:)
         integer, allocatable            :: globalToLocalElementID(:)
         
         integer                         :: bFaceOrder, numBFacePoints
         integer                         :: i, j, k, l, pNode, pElement
         integer                         :: fUnit, fileStat
         integer                         :: nodeIDs(NODES_PER_ELEMENT), nodeMap(NODES_PER_FACE)
         real(kind=RP)                   :: x(NDIM)
         integer                         :: faceFlags(FACES_PER_ELEMENT)
         CHARACTER(LEN=BC_STRING_LENGTH) :: names(FACES_PER_ELEMENT)
         TYPE(FacePatch), DIMENSION(6)   :: facePatches
         real(kind=RP)                   :: corners(NDIM,NODES_PER_ELEMENT)
!
!        ------------------
!        For curved patches
!        ------------------
!
         type(SurfInfo_t), allocatable                  :: SurfInfo(:)
         real(kind=RP)  , DIMENSION(:)    , ALLOCATABLE :: uNodes, vNodes
         real(kind=RP)  , DIMENSION(:,:,:), ALLOCATABLE :: values
!
!        ----------------
!        For flat patches
!        ----------------
!
         real(kind=RP)  , DIMENSION(2)     :: uNodesFlat = [-1.0_RP,1.0_RP]
         real(kind=RP)  , DIMENSION(2)     :: vNodesFlat = [-1.0_RP,1.0_RP]
         real(kind=RP)  , DIMENSION(3,2,2) :: valuesFlat
         character(len=LINE_LENGTH)        :: partitionName

         success               = .TRUE.
!
!        -----------------------
!        Read header information
!        -----------------------
!
         fUnit = UnusedUnit()
         OPEN( UNIT = fUnit, FILE = fileName, iostat = fileStat )
         IF ( fileStat /= 0 )     THEN
            PRINT *, "Error opening file: ", fileName
            success = .FALSE.
            RETURN 
         END IF
         
         READ(fUnit,*) numberOfAllNodes, numberOfAllElements, bFaceOrder

         self % nodeType = nodes
         self % no_of_elements = mpi_partition % no_of_elements
!
!        ----------------------------
!        Set up for face patches
!        Face patches are defined 
!        at chebyshev- lobatto points
!        ----------------------------
!
         numBFacePoints = bFaceOrder + 1
         allocate(uNodes(numBFacePoints))
         allocate(vNodes(numBFacePoints))
         allocate(values(3,numBFacePoints,numBFacePoints))
         
         DO i = 1, numBFacePoints
            uNodes(i) = -COS((i-1.0_RP)*PI/(numBFacePoints-1.0_RP)) 
            vNodes(i) = -COS((i-1.0_RP)*PI/(numBFacePoints-1.0_RP)) 
         END DO
         
         DO k = 1, 6 ! Most patches will be flat, so set up for self
            CALL facePatches(k) % construct(uNodesFlat,vNodesFlat) 
         END DO  
!
!        ---------------
!        Allocate memory
!        ---------------
!
         allocate( self % elements(mpi_partition % no_of_elements) )
         allocate( self % nodes(mpi_partition % no_of_nodes) )
         allocate( SurfInfo(mpi_partition % no_of_nodes) )
         allocate( globalToLocalNodeID(numberOfAllNodes) )
         allocate( globalToLocalElementID(numberOfAllElements) )
         self % no_of_elements = mpi_partition % no_of_elements
         globalToLocalNodeID = -1
         globalToLocalElementID = -1
!
!        ----------------------------------
!        Read nodes: Nodes have the format:
!        x y z
!        ----------------------------------
!
         pNode = 1
         DO j = 1, numberOfAllNodes 
            READ( fUnit, * ) x
            if ( pNode .gt. mpi_partition % no_of_nodes ) cycle
!
!           Construct only nodes that belong to the partition
!           -------------------------------------------------
            if ( j .eq. mpi_partition % nodeIDs(pNode) ) then
               CALL ConstructNode( self % nodes(pNode), x, j )
               globalToLocalNodeID(j) = pNode
               pNode = pNode + 1 
            end if
         END DO
!
!        -----------------------------------------
!        Read elements: Elements have the format:
!        node1ID node2ID node3ID node4ID ... node8ID
!        b1 b2 b3 b4 b5 b6
!           (=0 for straight side, 1 for curved)
!        If curved boundaries, then for each:
!           for j = 0 to bFaceOrder
!              x_j  y_j  z_j
!           next j
!        bname1 bname2 bname3 bname4 bname5 bname6
!        -----------------------------------------
!
         pElement = 1
         DO l = 1, numberOfAllElements 

            if ( pElement .gt. mpi_partition % no_of_elements ) then
!
!              Skip the element
!              ----------------
               READ( fUnit, * ) nodeIDs
               READ( fUnit, * ) faceFlags

               do k = 1, FACES_PER_ELEMENT
                  if ( faceFlags(k) .ne. 0 ) then
                     DO j = 1, numBFacePoints
                        DO i = 1, numBFacePoints
                           READ( fUnit, * )
                        END DO  
                     END DO
                  end if 
               end do
!
!              Get the boundary names to build zones (even they could be empty)
!              ----------------------------------------------------------------
               READ( fUnit, * ) names
               DO k = 1, 6
                  IF(TRIM(names(k)) /= emptyBCName) then
                     if ( all(trim(names(k)) .ne. zoneNameDictionary % allKeys()) ) then
                        call zoneNameDictionary % addValueForKey(trim(names(k)), trim(names(k)))
                     end if
                  end if
               END DO  

               cycle 

            else if ( l .ne. mpi_partition % elementIDs(pElement) ) then
!
!              Skip the element
!              ----------------
               READ( fUnit, * ) nodeIDs
               READ( fUnit, * ) faceFlags

               do k = 1, FACES_PER_ELEMENT
                  if ( faceFlags(k) .ne. 0 ) then
                     DO j = 1, numBFacePoints
                        DO i = 1, numBFacePoints
                           READ( fUnit, * )
                        END DO  
                     END DO
                  end if 
               end do
!
!              Get the boundary names to build zones (even they could be empty)
!              ----------------------------------------------------------------
               READ( fUnit, * ) names
               DO k = 1, 6
                  IF(TRIM(names(k)) /= emptyBCName) then
                     if ( all(trim(names(k)) .ne. zoneNameDictionary % allKeys()) ) then
                        call zoneNameDictionary % addValueForKey(trim(names(k)), trim(names(k)))
                     end if
                  end if
               END DO  

               cycle 
            end if
!
!           Otherwise, read the element
!           ---------------------------
            READ( fUnit, * ) nodeIDs
            READ( fUnit, * ) faceFlags
!
!           Translate nodeIDs to local nodeIDs
!           ----------------------------------
            nodeIDs = globalToLocalNodeID(nodeIDs)
!
!           -----------------------------------------------------------------------------
!           If the faceFlags are all zero, then self is a straight-sided
!           hex. In that case, set the corners of the hex8Map and use that in determining
!           the element geometry.
!           -----------------------------------------------------------------------------
!
            IF(MAXVAL(faceFlags) == 0)     THEN
               DO k = 1, NODES_PER_ELEMENT
                  corners(:,k) = self % nodes(nodeIDs(k)) % x
               END DO
               SurfInfo(pElement) % IsHex8 = .TRUE.
               SurfInfo(pElement) % corners = corners
               
            ELSE
!
!              --------------------------------------------------------------
!              Otherwise, we have to look at each of the faces of the element 
!              --------------------------------------------------------------
!
               DO k = 1, FACES_PER_ELEMENT 
                 IF ( faceFlags(k) == 0 )     THEN
!
!                    ----------
!                    Flat faces
!                    ----------
!
                     nodeMap           = localFaceNode(:,k)
                     valuesFlat(:,1,1) = self % nodes(nodeIDs(nodeMap(1))) % x
                     valuesFlat(:,2,1) = self % nodes(nodeIDs(nodeMap(2))) % x
                     valuesFlat(:,2,2) = self % nodes(nodeIDs(nodeMap(3))) % x
                     valuesFlat(:,1,2) = self % nodes(nodeIDs(nodeMap(4))) % x
                     
                     call SurfInfo(pElement) % facePatches(k) % construct(uNodesFlat, vNodesFlat, valuesFlat)
                     
                  ELSE
!
!                    -------------
!                    Curved faces 
!                    -------------
!
                     DO j = 1, numBFacePoints
                        DO i = 1, numBFacePoints
                           READ( fUnit, * ) values(:,i,j)
                        END DO  
                     END DO
                     
                     call SurfInfo(pElement) % facePatches(k) % construct(uNodes, vNodes, values)
                     
                  END IF
               END DO
               
            END IF 
!
!           -------------------------
!           Now construct the element
!           -------------------------
!
            call self % elements(pElement) % Construct (Nx(l), Ny(l), Nz(l), nodeIDs , pElement, l)
            
            READ( fUnit, * ) names
            CALL SetElementBoundaryNames( self % elements(pElement), names )
!
!           Add the boundary face to the zone Name dictionary
!           -------------------------------------------------            
            DO k = 1, 6
               IF(TRIM(names(k)) /= emptyBCName) then
                  if ( all(trim(names(k)) .ne. zoneNameDictionary % allKeys()) ) then
                     call zoneNameDictionary % addValueForKey(trim(names(k)), trim(names(k)))
                  end if
               end if
            END DO  
!
!           Count the partition element
!           ---------------------------            
            globalToLocalElementID(l) = pElement
            pElement = pElement + 1
            
         END DO      ! l = 1, numberOfAllElements
!
!        ---------------------------
!        Construct the element faces
!        ---------------------------
!
         numberOfFaces        = GetOriginalNumberOfFaces(self)
         self % numberOfFaces = numberOfFaces
         allocate( self % faces(self % numberOfFaces) )
         CALL ConstructFaces( self, success )
!
!        --------------------------------
!        Get actual mesh element face IDs
!        --------------------------------
!
         CALL getElementsFaceIDs(self)
!
!        --------------
!        Cast MPI faces
!        --------------
!
         call ConstructMPIFaces
         call self % UpdateFacesWithPartition(mpi_partition, &
                                              numberOfAllElements, &
                                              globalToLocalElementID)
!
!        -------------------------
!        Build the different zones
!        -------------------------
!
         call self % ConstructZones()
!
!        ---------------------------
!        Construct periodic faces
!        ---------------------------
!
         CALL ConstructPeriodicFaces( self )
!
!        ---------------------------
!        Delete periodic- faces
!        ---------------------------
!
         CALL DeletePeriodicMinusFaces( self )
!
!        ---------------------------
!        Assign faces ID to elements
!        ---------------------------
!
         CALL getElementsFaceIDs(self)
!
!        ---------------------
!        Define boundary faces
!        ---------------------
!
         call self % DefineAsBoundaryFaces()
!
!        -----------------------------------
!        Check if this is a 2D extruded mesh
!        -----------------------------------
!
         call self % CheckIfMeshIs2D()
!
!        -------------------------------
!        Set the mesh as 2D if requested
!        -------------------------------
!
         if ( dir2D .ne. 0 ) then
            call SetMappingsToCrossProduct
            call self % CorrectOrderFor2DMesh(dir2D)
         end if
!
!        ------------------------------
!        Set the element connectivities
!        ------------------------------
!
         call self % SetConnectivitiesAndLinkFaces(nodes)
!
!        ---------------------------------------
!        Construct elements' and faces' geometry
!        ---------------------------------------
!
         call self % ConstructGeometry(SurfInfo)

         CLOSE( fUnit )
!
!        ---------
!        Finish up
!        ---------
!
         if (.not. self % child) CALL self % DescribePartition( trim(fileName) )
!
!        --------------------
!        Prepare mesh for I/O
!        --------------------
!
         call self % PrepareForIO
         if (.not. self % child) call self % Export( trim(fileName) )

         deallocate(globalToLocalNodeID)
         deallocate(globalToLocalElementID)

      END SUBROUTINE ConstructMeshPartition_FromSpecMeshFile_

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      function NumOfElems_SpecMesh( fileName ) result(nelem)
         implicit none
         !----------------------------------
         CHARACTER(LEN=*), intent(in) :: fileName
         integer                      :: nelem
         !----------------------------------
         integer :: k, fUnit
         !----------------------------------
         
         OPEN(newunit = fUnit, FILE = trim(fileName) )  
            READ(fUnit,*) k, nelem, k                    ! Here k is used as default reader since this variables are not important now
         CLOSE(fUnit)
      
      end function NumOfElems_SpecMesh
end module Read_SpecMesh
