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
      MODULE HexMeshClass
      USE MeshTypes
      USE NodeClass
      USE ElementClass
      USE FaceClass
      USE TransfiniteMapClass
      use SharedBCModule
      IMPLICIT NONE
!
!     ---------------
!     Mesh definition
!     ---------------
!
      TYPE HexMesh
         INTEGER                                  :: numberOfFaces
         TYPE(Node)   , DIMENSION(:), ALLOCATABLE :: nodes
         TYPE(Face)   , DIMENSION(:), ALLOCATABLE :: faces, dummy_faces
         TYPE(Element), DIMENSION(:), ALLOCATABLE :: elements
!
!        ========         
         CONTAINS
!        ========         
!
         PROCEDURE :: constructFromFile => ConstructMesh_FromFile_
         PROCEDURE :: destruct          => DestructMesh
         
      END TYPE HexMesh
!
!     ========
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructMesh_FromFile_( self, fileName, spA, success )
         IMPLICIT NONE
!
!        ---------------
!        Input variables
!        ---------------
!
         CLASS(HexMesh)     :: self
         TYPE(NodalStorage) :: spA
         CHARACTER(LEN=*)   :: fileName
         LOGICAL            :: success
!
!        ---------------
!        Local variables
!        ---------------
!
         TYPE(TransfiniteHexMap), POINTER :: hexMap, hex8Map, genHexMap
         
         INTEGER                         :: N
         INTEGER                         :: numberOfElements
         INTEGER                         :: numberOfNodes
         INTEGER                         :: numberOfBoundaryFaces
         INTEGER                         :: numberOfFaces
         
         INTEGER                         :: bFaceOrder, numBFacePoints
         INTEGER                         :: i, j, k, l
         INTEGER                         :: fUnit, fileStat
         INTEGER                         :: nodeIDs(8), nodeMap(4)
         REAL(KIND=RP)                   :: x(3)
         INTEGER                         :: faceFlags(6)
         CHARACTER(LEN=BC_STRING_LENGTH) :: names(6)
         TYPE(FacePatch), DIMENSION(6)   :: facePatches
         REAL(KIND=RP)                   :: corners(3,8)
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
!
!        ---------------
!        For Periodic BC
!        ---------------
!
         REAL(KIND=RP) :: x1(3), x2(3)
         LOGICAL       :: success_periodical(4), flag_periodical(4)
         INTEGER       :: coord 
          
         numberOfBoundaryFaces = 0
         N                     = spa % N
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
!        b1 b2 b3 b4 ... b8
!           (=0 for straight side, 1 for curved)
!        If curved boundaries, then for each:
!           for j = 0 to bFaceOrder
!              x_j  y_j z_j
!           next j
!        bname1 bname2 bname3 bname4 ... bname8
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
               DO k = 1, 8
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
               DO k = 1, 6 
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
                        
                     IF(facePatches(k) % noOfKnots(1) == 2)     THEN
                        CALL facePatches(k) % destruct()
                        CALL facePatches(k) % construct(uNodes, vNodes, values)
                     ELSE
                       CALL facePatches(k) % setFacePoints(points = values)
                     END IF 
                     
                  END IF
               END DO
               CALL genHexMap % destruct()
               CALL genHexMap % constructWithFaces(facePatches)
               
               hexMap => genHexMap
            END IF 
!
!           -------------------------
!           Now construct the element
!           -------------------------
!
            
            CALL ConstructElementGeometry( self % elements(l), spA, nodeIDs, hexMap )
            
            READ( fUnit, * ) names
            CALL SetElementBoundaryNames( self % elements(l), names )
            
            DO k = 1, 6
               IF(TRIM(names(k)) /= emptyBCName) numberOfBoundaryFaces = numberOfBoundaryFaces + 1
            END DO  
         END DO
!
!        ---------------------------
!        Construct the element faces
!        ---------------------------
!
         numberOfFaces        = (6*numberOfElements + numberOfBoundaryFaces)/2
         self % numberOfFaces = numberOfFaces
         ALLOCATE( self % faces(self % numberOfFaces) )
         CALL ConstructFaces( self, success )

         DO i = 1, self%numberOfFaces
            PRINT*, i, self % faces(i) % elementSide(1),  "facei"
            PRINT*, self % faces(i) % boundaryName
            PRINT*, TRIM(bcTypeDictionary % stringValueForKey(key             = self%faces(i)%boundaryName, &
                                                         requestedLength = BC_STRING_LENGTH))
            IF (TRIM(bcTypeDictionary % stringValueForKey(key             = self%faces(i)%boundaryName, &
                                                         requestedLength = BC_STRING_LENGTH)) == "periodic+") THEN  
               DO j = 1, self%numberOfFaces
                  PRINT*, j, self % faces(j) % elementSide(1), "facej"     
                  PRINT*, self % faces(j) % boundaryName    
                  PRINT*, trim(bcTypeDictionary % stringValueForKey(key             = self%faces(j)%boundaryName, &
                                                         requestedLength = BC_STRING_LENGTH))   
                  IF ((trim(bcTypeDictionary % stringValueForKey(key             = self%faces(j)%boundaryName, &
                                                         requestedLength = BC_STRING_LENGTH)) == "periodic-")) THEN  
                     coord = 0
                     success_periodical(:) = .FALSE.
                     flag_periodical(:)    = .FALSE.
                     DO k = 1, 4
                        PRINT*, "k", k
                        x1 = self%nodes(self%faces(i)%nodeIDs(k))%x
                        DO l = 1, 4
                           IF (.NOT.flag_periodical(l)) THEN 
                              PRINT*, "l", l
                              x2 = self%nodes(self%faces(j)%nodeIDs(l))%x   
                              PRINT*, "x1", x1
                              PRINT*, "x2", x2                    
                              CALL CompareTwoNodes(x1, x2, success_periodical(k), coord)
                              PRINT*, success_periodical(k), coord
                              PRINT*, "----------------------------------------"
                              IF (success_periodical(k)) THEN 
                              !   PRINT*, k,l, coord
                              !   PRINT*, x1
                              !   PRINT*, x2
                              !   PRINT*, success_periodical
                              !   PRINT*, "-------------------------------"
                                 flag_periodical(l) = .TRUE. 
                                 EXIT
                              ENDIF  
                           ENDIF 
                        !Estoy hay que hacer que se comparen todos los nodos dos a dos
                        !En función de como coincidan, hay que organizar la rotacion de las caras
                        ENDDO 
                        IF (.NOT.success_periodical(k)) EXIT  
                     ENDDO  
                        PRINT*,         success_periodical                 
                        IF ((success_periodical(1)).AND.(success_periodical(2))&
                        .AND.(success_periodical(3)).AND.(success_periodical(4))) THEN
                        PRINT*, "I'm IN",        success_periodical                 
                        
                           !DO l = 1,4
                           !   self%nodes(self%faces(j)%nodeIDs(l))%x(coord) = self%nodes(self%faces(i)%nodeIDs(l))%x(coord)
                           !ENDDO 
!                           PRINT*, "x1 ", x1
!                           PRINT*, "x2 ", x2
!                           PRINT*, "coord", coord 
!                           PRINT*, "i,j ", i,j
!                           PRINT*, "____________________________________"

                           self % faces(i) % boundaryName = ""
                           self % faces(i) % elementIDs(2)  = self % faces(j) % elementIDs(1)
                           self % faces(i) % elementSide(2) = self % faces(j) % elementSide(1) 
                           self % faces(i) % FaceType       = HMESH_INTERIOR
                           self % faces(i) % rotation       = 0!faceRotation(masterNodeIDs = self % faces(i) % nodeIDs, &
                                                              !           slaveNodeIDs  = self % faces(i) % nodeIDs)      
                                                                                     
                        ENDIF    

                  ENDIF 
               ENDDO
            ENDIF 
         ENDDO     
         j = 0
         ALLOCATE( self % dummy_faces(self % numberOfFaces) )
         DO i = 1, self%numberOfFaces 
            IF (trim(bcTypeDictionary % stringValueForKey(key             = self%faces(i)%boundaryName, &
                                                         requestedLength = BC_STRING_LENGTH)) /= "periodic-") THEN 
               j = j+1
               self%dummy_faces(j) = self%faces(i)
            ENDIF 
         ENDDO 
         
         DEALLOCATE(self%faces)
         ALLOCATE(self%faces(j))
         self%numberOfFaces = j
         DO i = 1, self%numberOfFaces
            self%faces(i) = self%dummy_faces(i)
         ENDDO 
         
               !Una vez encontrada, encuentro las coordenadas de sus nodos
               !Busco otra cara (nuevo bucle) que sea periodica Y cuyas coordenadas de nodos (salvo una) coincidan
            
               !Cuando encuentro una pareja, cambio la condición de contorno de ambas, las pongo como condicion de contorno anterior y le asigno los valores de la encontrada
!            self % faces(self % numberOfFaces) % boundaryName = "---"
!            self % faces(faceID) % elementIDs(2)  = eID
!            self % faces(faceID) % elementSide(2) = faceNumber
!            self % faces(faceID) % FaceType       = HMESH_INTERIOR
!            self % faces(faceID) % rotation       = faceRotation(masterNodeIDs = self % faces(faceID) % nodeIDs, &
!                                                                       slaveNodeIDs  = faceNodeIDs)
            !Puedo hacer una comprobación de que no queda ninguna "Periodica +"
         
         !Ahora habría que volver a construir las caras, sabiendo cual es el número correcto de las mismas
         !Estos es facil de saber a priori, descontando las periodic - de la cuenta incial
         !Se van copiando todas de la estructura inicial a la nueva, dejando fuera a las que tengan como condicion de 
         !contorno periodic -
            
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
         
      END SUBROUTINE ConstructMesh_FromFile_
      
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
!-------------------- 
! External variables 
!-------------------- 
!  
      REAL(KIND=RP) :: x1(3)
      REAL(KIND=RP) :: x2(3)
      LOGICAL       :: success
      INTEGER       :: coord 
! 
!-------------------- 
! Local variables 
!-------------------- 
! 
      INTEGER :: i
      INTEGER :: counter    
      
      counter = 0
      
      IF (coord == 0) THEN

         DO i = 1,3
            IF (x1(i) == x2(i)) THEN 
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
               IF (x1(i) == x2(i)) THEN 
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
          
             
      END SUBROUTINE  
! 
!//////////////////////////////////////////////////////////////////////// 
! 
      !^ 
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
            CALL DestructFace( self % faces(j) )
         END DO
         DEALLOCATE( self % faces )
         
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
         CALL PrintFace( self % faces(k))
      END DO
      
      END SUBROUTINE PrintMesh
!
!//////////////////////////////////////////////////////////////////////// 
! 
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
         
         TYPE(FTMultiIndexTable), POINTER  :: table
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
                                                                       slaveNodeIDs  = faceNodeIDs)
               ELSE 
!
!                 ------------------
!                 Construct new face
!                 ------------------
!
                  self % numberOfFaces = self % numberOfFaces + 1
                  
                  IF(self % numberOfFaces > SIZE(self % faces))     THEN
          
                     CALL table % release()
                     DEALLOCATE(table)
                     PRINT *, "Too many faces for # of elements:", self % numberOfFaces, " vs ", SIZE(self % faces)
                     success = .FALSE.
                     RETURN  
                  END IF 
                  
                  CALL ConstructFace( self % faces(self % numberOfFaces),   &
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
                  CALL v % release()
               END IF 
            END DO 
              
         END DO  
         
         CALL table % release()
         DEALLOCATE(table)
         
      END SUBROUTINE ConstructFaces
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
      INTEGER FUNCTION faceRotation( masterNodeIDs, slaveNodeIDs)
         IMPLICIT NONE 
         INTEGER, DIMENSION(4) :: masterNodeIDs, slaveNodeIDs
         
         INTEGER :: j
         
         DO j = 1, 4
            IF(masterNodeIDs(1) == slaveNodeIDs(j)) EXIT 
         END DO  
         faceRotation = j - 1
          
      END FUNCTION faceRotation
      
!     ==========
      END MODULE HexMeshClass
!     ==========
!
      