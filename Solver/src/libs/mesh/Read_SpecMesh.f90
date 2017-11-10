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
      IMPLICIT NONE
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
         TYPE(TransfiniteHexMap), POINTER :: hexMap, hex8Map, genHexMap
         
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
         
         real(kind=RP), allocatable, dimension(:)     :: xiCL,etaCL,zetaCL
         real(kind=RP), allocatable, dimension(:,:,:) :: face1CL,face2CL,face3CL,face4CL,face5CL,face6CL
          
         numberOfBoundaryFaces = 0
         corners               = 0.0_RP
         success               = .TRUE.
         
         ALLOCATE(hex8Map)
         CALL hex8Map % constructWithCorners(corners)
         ALLOCATE(genHexMap)
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
            uNodes(i) = -COS((i-1)*PI/(numBFacePoints-1)) 
            vNodes(i) = -COS((i-1)*PI/(numBFacePoints-1)) 
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
               CALL hex8Map % setCorners(corners)
               hexMap => hex8Map
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
                     
                     IF(facePatches(k) % noOfKnots(1) /= 2)     THEN
                        CALL facePatches(k) % destruct()
                        CALL facePatches(k) % construct(uNodesFlat, vNodesFlat, valuesFlat)
                     ELSE
                        CALL facePatches(k) % setFacePoints(points = valuesFlat)
                     END IF 
                     
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
                        
                     IF(facePatches(k) % noOfKnots(1) /= numBFacePoints)     THEN      ! TODO This could be problematic with anisotropy
                        CALL facePatches(k) % destruct()
                        CALL facePatches(k) % construct(uNodes, vNodes, values)
                     ELSE
                       CALL facePatches(k) % setFacePoints(points = values)
                     END IF 
                     
                  END IF
               END DO
               
!
!              Impose subparametric, or at most isoparametric, representation of boundaries
!              To do so, we interpolate the boundary points to (Nx+1, Ny+1, Nz+1) 
!              Chebyshev-Lobatto nodes and reconstruct the patches there. 
!              TODO: This can cause problems with p-adaptation for inner curved boundaries
!              ----------------------------------------------------------------------------
            
               ! Allocation
               
               allocate(xiCL(Nx(l)+1),etaCL(Ny(l)+1),zetaCL(Nz(l)+1))
               allocate(face1CL(1:3,Nx(l)+1,Nz(l)+1))
               allocate(face2CL(1:3,Nx(l)+1,Nz(l)+1))
               allocate(face3CL(1:3,Nx(l)+1,Ny(l)+1))
               allocate(face4CL(1:3,Ny(l)+1,Nz(l)+1))
               allocate(face5CL(1:3,Nx(l)+1,Ny(l)+1))
               allocate(face6CL(1:3,Ny(l)+1,Nz(l)+1))
            
               ! Construct the interpolants based on Chebyshev-Lobatto points
               
               xiCL   = (/ ( -cos((i-1)*PI/Nx(l)),i=1, Nx(l)+1) /)
               etaCL  = (/ ( -cos((i-1)*PI/Ny(l)),i=1, Ny(l)+1) /)
               zetaCL = (/ ( -cos((i-1)*PI/Nz(l)),i=1, Nz(l)+1) /)
               
               ! Project the faces to new boundary representations
               
               call ProjectFaceToNewPoints(facePatches(1), Nx(l), xiCL , Nz(l), zetaCL, face1CL)
               call ProjectFaceToNewPoints(facePatches(2), Nx(l), xiCL , Nz(l), zetaCL, face2CL)
               call ProjectFaceToNewPoints(facePatches(3), Nx(l), xiCL , Ny(l), etaCL , face3CL)
               call ProjectFaceToNewPoints(facePatches(4), Ny(l), etaCL, Nz(l), zetaCL, face4CL)
               call ProjectFaceToNewPoints(facePatches(5), Nx(l), xiCL , Ny(l), etaCL , face5CL)
               call ProjectFaceToNewPoints(facePatches(6), Ny(l), etaCL, Nz(l), zetaCL, face6CL)
               
               ! Destruct face patches
               
               call facePatches(1) % Destruct()
               call facePatches(2) % Destruct()
               call facePatches(3) % Destruct()
               call facePatches(4) % Destruct()
               call facePatches(5) % Destruct()
               call facePatches(6) % Destruct()
               
               ! Construct the new Chebyshev-Lobatto face patches
               
               call facePatches(1) % Construct( xiCL , zetaCL, face1CL )
               call facePatches(2) % Construct( xiCL , zetaCL, face2CL )
               call facePatches(3) % Construct( xiCL , etaCL , face3CL )
               call facePatches(4) % Construct( etaCL, zetaCL, face4CL )
               call facePatches(5) % Construct( xiCL , etaCL , face5CL )
               call facePatches(6) % Construct( etaCL, zetaCL, face6CL )
               
               deallocate(xiCL,etaCL,zetaCL)
               deallocate(face1CL,face2CL,face3CL,face4CL,face5CL,face6CL)
         
!
!              Construct the transfinite map with the boundary representations
!              ---------------------------------------------------------------
               
               CALL genHexMap % destruct()
               CALL genHexMap % constructWithFaces(facePatches)
               
               hexMap => genHexMap
            END IF 
!
!           -------------------------
!           Now construct the element
!           -------------------------
!
            CALL ConstructElementGeometry( self % elements(l), spA(Nx(l)), spA(Ny(l)), spA(Nz(l)), nodeIDs, hexMap , l)
            
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
            
!
!           ------------------------------------
!           Construct the element connectivities
!           ------------------------------------
!
            DO k = 1, 6
               IF (TRIM(names(k)) == emptyBCName ) THEN
                  self%elements(l)%NumberOfConnections(k) = 1
                  CALL self%elements(l)%Connection(k)%construct (1)  ! Just conforming elements
               ELSE
                  self%elements(l)%NumberOfConnections(k) = 0
               ENDIF
            ENDDO
            
            
         END DO      ! l = 1, numberOfElement
        
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
!        ------------------------------
!        Set the element connectivities
!        ------------------------------
!
         call self % SetConnectivities
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
            
         CLOSE( fUnit )
!
!        ---------
!        Finish up
!        ---------
!
         CALL hex8Map % destruct()
         DEALLOCATE(hex8Map)
         CALL genHexMap % destruct()
         DEALLOCATE(genHexMap)

         CALL self % Describe( trim(fileName) )
         
         self % Ns = Nx

         call self % Export( trim(fileName) )
         
      END SUBROUTINE ConstructMesh_FromSpecMeshFile_
      
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
