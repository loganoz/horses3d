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
      use HexMeshClass
      use SMConstants
      USE TransfiniteMapClass
      use FacePatchClass
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
      SUBROUTINE ConstructMesh_FromSpecMeshFile_( self, fileName, nodes, spA, Nx, Ny, Nz, success )
         USE Physics
         IMPLICIT NONE
!
!        ---------------
!        Input variables
!        ---------------
!
         CLASS(HexMesh)     :: self
         integer            :: nodes
         TYPE(NodalStorage) :: spA(0:)
         CHARACTER(LEN=*)   :: fileName
         INTEGER            :: Nx(:), Ny(:), Nz(:)     !<  Polynomial orders for all the elements
         LOGICAL            :: success
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                         :: numberOfElements
         INTEGER                         :: numberOfNodes
         INTEGER                         :: numberOfBoundaryFaces
         INTEGER                         :: numberOfFaces
         
         INTEGER                         :: bFaceOrder, numBFacePoints
         INTEGER                         :: i, j, k, l
         INTEGER                         :: fUnit, fileStat
         INTEGER                         :: nodeIDs(NODES_PER_ELEMENT), nodeMap(NODES_PER_FACE)
         REAL(KIND=RP)                   :: x(NDIM)
         INTEGER                         :: faceFlags(FACES_PER_ELEMENT)
         CHARACTER(LEN=BC_STRING_LENGTH) :: names(FACES_PER_ELEMENT)
         TYPE(FacePatch), DIMENSION(6)   :: facePatches
         REAL(KIND=RP)                   :: corners(NDIM,NODES_PER_ELEMENT)
         INTEGER, EXTERNAL               :: UnusedUnit
!
!        ------------------
!        For curved patches
!        ------------------
!
         type(SurfInfo_t), allocatable                  :: SurfInfo(:)
         REAL(KIND=RP)  , DIMENSION(:)    , ALLOCATABLE :: uNodes, vNodes
         REAL(KIND=RP)  , DIMENSION(:,:,:), ALLOCATABLE :: values
!
!        ----------------
!        For flat patches
!        ----------------
!
         REAL(KIND=RP)  , DIMENSION(2)     :: uNodesFlat = [-1.0_RP,1.0_RP]
         REAL(KIND=RP)  , DIMENSION(2)     :: vNodesFlat = [-1.0_RP,1.0_RP]
         REAL(KIND=RP)  , DIMENSION(3,2,2) :: valuesFlat
          
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
!        ---------------
!        Allocate memory
!        ---------------
!
         ALLOCATE( self % elements(numberOfelements) )
         ALLOCATE( self % nodes(numberOfNodes) )
         allocate( SurfInfo(numberOfelements) )
!
!        ----------------------------
!        Set up for face patches
!        Face patches are defined 
!        at chebyshev- lobatto points
!        ----------------------------
!
         numBFacePoints = bFaceOrder + 1
         ALLOCATE(uNodes(numBFacePoints))
         ALLOCATE(vNodes(numBFacePoints))
         ALLOCATE(values(3,numBFacePoints,numBFacePoints))
         
         DO i = 1, numBFacePoints
            uNodes(i) = -COS((i-1.0_RP)*PI/(numBFacePoints-1.0_RP)) 
            vNodes(i) = -COS((i-1.0_RP)*PI/(numBFacePoints-1.0_RP)) 
         END DO
         
         DO k = 1, 6 ! Most patches will be flat, so set up for self
            CALL facePatches(k) % construct(uNodesFlat,vNodesFlat) 
         END DO  
!
!        ----------------------------------
!        Read nodes: Nodes have the format:
!        x y z
!        ----------------------------------
!
         DO j = 1, numberOfNodes 
            READ( fUnit, * ) x
            CALL ConstructNode( self % nodes(j), x )
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
            call self % elements(l) % Construct (spA(Nx(l)), spA(Ny(l)), spA(Nz(l)), nodeIDs , l)
            
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
         ALLOCATE( self % faces(self % numberOfFaces) )
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
!        ------------------------------
!        Set the element connectivities
!        ------------------------------
!
         call self % SetConnectivities(spA,nodes)
!
!        ---------------------------------------
!        Construct elements' and faces' geometry
!        ---------------------------------------
!
         call self % ConstructGeometry(spA,SurfInfo)
            
         CLOSE( fUnit )
!
!        ---------
!        Finish up
!        ---------
!

         CALL self % Describe( trim(fileName) )
         
         self % Ns = Nx

         call self % Export( trim(fileName) )
         
      END SUBROUTINE ConstructMesh_FromSpecMeshFile_
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
