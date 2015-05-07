!
!////////////////////////////////////////////////////////////////////////
!
!      QuadMesh.f95
!      Created: 2007-03-22 17:05:00 -0400 
!      By: David Kopriva  
!
!      Implements Algorithms:
!         Algorithm 126: QuadMesh
!         Algorithm 127: ConstructMesh_FromFile_ (Construct)
!         Algorithm 148: ConstructEdges (ConstructMeshEdges)
!
!////////////////////////////////////////////////////////////////////////
!
      MODULE QuadMeshClass
      USE MeshTypes
      USE NodeClass
      USE ElementClass
      USE EdgeClass
      USE LinkedListClass
      USE SparseMatrixClass
      USE Nodal2DStorageClass
      USE TransfiniteMapClass
      IMPLICIT NONE
!
!     ---------------
!     Mesh definition
!     ---------------
!
      TYPE QuadMesh
         INTEGER                                  :: numberOfEdges
         TYPE(Edge)   , DIMENSION(:), ALLOCATABLE :: edges
         TYPE(node)   , DIMENSION(:), ALLOCATABLE :: nodes
         TYPE(Element), DIMENSION(:), ALLOCATABLE :: elements
      END TYPE QuadMesh
      
      INTEGER, DIMENSION(2,4) :: cornerMap
      INTEGER, DIMENSION(4)   :: sideMap
!
!     ----------
!     Interfaces
!     ----------
!
      INTERFACE Destruct
         MODULE PROCEDURE DestructMesh
      END INTERFACE Destruct
!
!     ========
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructMesh_FromFile_( this, spA, fileName, N )
         IMPLICIT NONE
!
!        ---------------
!        Input variables
!        ---------------
!
         INTEGER              :: N
         TYPE(QuadMesh)       :: this
         TYPE(Nodal2DStorage) :: spA
         CHARACTER(LEN=*)     :: fileName
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                         :: numberOfElements, numberOfNodes, bCurveOrder
         INTEGER                         :: i, j, k
         INTEGER                         :: fUnit = 4, stat
         INTEGER                         :: nodeIDs(4)
         REAL(KIND=RP)                   :: x(2), xStart(2), xEnd(2)
         INTEGER                         :: nEdges
         TYPE(TransfiniteQuadMap)        :: quadMap
         INTEGER                         :: curveFlags(4)
         CHARACTER(LEN=BC_STRING_LENGTH) :: names(4)
         
         TYPE(CurveInterpolant)       , POINTER     :: boundaryCurves(:)
         REAL(KIND=RP), DIMENSION(:)  , ALLOCATABLE :: nodes
         REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE :: values
!
!        -----------------------
!        Read header information
!        -----------------------
!
         OPEN( UNIT = fUnit, FILE = fileName, iostat = stat )
         IF ( STAT /= 0 )     THEN
            PRINT *, "Error opening file: ", fileName
            STOP
         END IF
         
         READ(fUnit,*) numberOfNodes, numberOfElements, bCurveOrder
!
!        ---------------
!        Allocate memory
!        ---------------
!
         ALLOCATE( this%elements(numberOfelements) )
         ALLOCATE( this%nodes(numberOfNodes) )
!
!        -------------------------------
!        Compute memory needed for edges
!        -------------------------------
!
         nEdges = numberOfelements + 4*numberOfNodes/3 - 1
         ALLOCATE( this%edges(nEdges) )
         this%numberOfEdges = 0
         
         ALLOCATE( boundaryCurves(4) )
         ALLOCATE( nodes(0:bCurveOrder) )
         ALLOCATE( values(0:bCurveOrder,2) )
!
!        -----------
!        Set globals
!        -----------
!
         cornerMap(:,1) = (/ 0,0 /)
         cornerMap(:,2) = (/ N,0 /)
         cornerMap(:,3) = (/ N,N /)
         cornerMap(:,4) = (/ 0,N /)
         
         sideMap(1) = 0
         sideMap(2) = N
         sideMap(3) = N
         sideMap(4) = 0
!
!        ----------------------
!        Set up boundary curves
!        and mapping
!        ----------------------
!
         DO j = 0, bCurveOrder 
            nodes(j) = -COS(j*PI/bCurveOrder)
         END DO
         values = 0.0_RP
         DO j = 1, 4 
            CALL Construct( boundaryCurves(j), bCurveOrder, nodes, values )
         END DO
         quadMap = NewTransfiniteQuadMap( boundaryCurves, MAP_DOESNT_OWN_CURVES )
!
!        ---------------------------------
!        Read nodes: Nodes have the format
!         x y
!        ---------------------------------
!
         DO j = 1, numberOfNodes 
            READ( fUnit, * ) x
            CALL ConstructNode( this%nodes(j), x )
         END DO
!
!        -----------------------------------------
!        Read elements: Elements have the format
!        node1 node2 node3 node4
!        b1 b2 b3 b4
!        (=0 for straight side, 1 for curved)
!        if curved boundaries, then for each:
!        for j = 0 to bCurveOrder
!           x_j  y_j
!        next j
!        bname1 bname2 bname3 bname4
!        -----------------------------------------
!
         DO j = 1, numberOfElements 
            READ( fUnit, * ) nodeIDs
            READ( fUnit, * ) curveFlags
            DO k = 1, 4 
               IF ( curveFlags(k) == 0 )     THEN
                  xStart = this%nodes(nodeIDs(edgeMap(1,k)))%x
                  xEnd   = this%nodes(nodeIDs(edgeMap(2,k)))%x
                  CALL ComputeLinearInterpolant( xStart, xEnd, bCurveOrder, nodes, values )
               ELSE
                  DO i = 0, bCurveOrder 
                     READ( fUnit, * ) values(i,:)
                  END DO
               END IF
               CALL SetValues( quadMap%boundaryCurves(k), values )
            END DO
            CALL ConstructElement( this%elements(j), spA, nodeIDs, quadMap )
            READ( fUnit, * ) names
            CALL SetElementBoundaryNames( this%elements(j), names )
         END DO
         
         CLOSE( fUnit )
!
!        ---------
!        Finish up
!        ---------
!
         CALL ConstructEdges( this )
         CALL FinishEdges( this, N )
         CALL MakeNodeToElementConnections( this )
         CALL SetBoundaryTypes( this )
         
         CALL Destruct( quadMap )
         DO k = 1, 4 
            CALL Destruct( boundaryCurves(k) )
         END DO
         DEALLOCATE( nodes, values, boundaryCurves )
         
      END SUBROUTINE ConstructMesh_FromFile_
!
!     -----------
!     Destructors
!     -----------
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructMesh( this )
         IMPLICIT NONE 
         TYPE(QuadMesh) :: this
         INTEGER        :: j
!
!        -----
!        Nodes
!        -----
!
         DO j = 1, SIZE(this%nodes) 
            CALL Destruct(this%nodes(j)%adjElements)
         END DO
         DEALLOCATE( this%nodes )
!
!        -----
!        Edges
!        -----
!
         DO j = 1, SIZE(this%edges) 
            CALL DestructEdge(this%edges(j))
         END DO
         DEALLOCATE( this%edges )
!
!        --------
!        Elements
!        --------
!
         DO j = 1, SIZE(this%elements) 
            CALL DestructElement( this%elements(j) )
         END DO
         DEALLOCATE( this%elements )

         
      END SUBROUTINE DestructMesh
!
!     -------------
!     Other methods
!     -------------
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE MakeNodeToElementConnections( this ) 
         IMPLICIT NONE
         TYPE(QuadMesh)     :: this
         INTEGER            :: e, n
         TYPE(ListData)     :: d
         
         DO e = 1, SIZE(this%elements) 
           DO n = 1,4
               CALL Construct( d, e, cornerMap(1,n), cornerMap(2,n) )
               CALL Add( this%nodes(this%elements(e)%nodeIDs(n))%adjElements, d  )
            END DO 
         END DO
      END SUBROUTINE MakeNodeToElementConnections
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructEdges( this ) 
         IMPLICIT NONE 
         TYPE(QuadMesh) :: this
         INTEGER        :: j, k, endNodes(2), s1, e1, n1, edgeID
         INTEGER        :: nodeIDs(4)
!
!        --------------
!        For hash table
!        --------------
!
         TYPE(SparseMatrix) :: table
         INTEGER            :: keys(2)
         INTEGER, EXTERNAL  :: Hash1, Hash2
         TYPE(ListData)     :: d
         
         CALL Construct( table, SIZE(this%nodes) )
         
         DO j = 1, SIZE(this%elements)
            nodeIDs = this%elements(j)%nodeIDs
            DO k = 1,4
               endNodes = (/ nodeIDs(edgeMap(1,k)), nodeIDs(edgeMap(2,k)) /)
               keys(1) = Hash1(endNodes)
               keys(2) = Hash2(endNodes)
               
               IF ( ContainsKeys( table, keys(1), keys(2) ) )     THEN
                  CALL DataForKeys( table, keys(1), keys(2), d )
                  edgeID = d%id
!
!                 ---------------------
!                 Original element info
!                 ---------------------
!
                  s1 = this%edges(edgeID)%elementSide(1)
                  e1 = this%edges(edgeID)%elementIDs(1)
                  n1 = this%elements(e1)%nodeIDs( (edgeMap(1,s1)) )
!
!                 ------------------------
!                 Set current element info
!                 ------------------------
!
                  this%edges(edgeID)%elementIDs(2)  = j
                  IF ( endNodes(1) == n1 )     THEN
                     this%edges(edgeID)%elementSide(2) = k
                  ELSE
                     this%edges(edgeID)%elementSide(2) = -k
                  END IF
               ELSE
                  this%numberOfEdges = this%numberOfEdges + 1
                  edgeID = this%numberOfEdges
                  CALL ConstructEdge( this%edges(edgeID), endNodes, j, k )
                  CALL Construct( d, edgeID, 0, 0 )
                  CALL AddDataForKeys( table, keys(1), keys(2), d )
               END IF
               
            END DO
         END DO
!
!        -------
!        Cleanup
!        -------
!
         CALL Destruct( table ) ! not needed anymore
         
      END SUBROUTINE ConstructEdges
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SetBoundaryTypes( this ) 
         IMPLICIT NONE 
         TYPE(QuadMesh)                  :: this
         INTEGER                         :: j, elementID, elementSide, k
         CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryName
         
         DO j = 1, this%numberOfEdges 
            IF( this%edges(j)%elementIDs(2) == QMESH_NONE )     THEN
!
!              ---------------------------------
!              Set the edge boundary information
!              ---------------------------------
!
               this%edges(j)%edgeType = QMESH_BOUNDARY
               elementID   = this%edges(j)%elementIDs(1)
               elementSide = this%edges(j)%elementSide(1)
               boundaryName               = this%elements(elementID)%boundaryName(elementSide)
               this%edges(j)%boundaryName = boundaryName
!
!              ---------------------------------
!              Set the node boundary information
!              ---------------------------------
!
               DO k = 1, 2 
                  this%nodes(this%edges(j)%nodeIds(k))%nodeType = QMESH_BOUNDARY
                  
                  IF ( this%nodes(this%edges(j)%nodeIds(k))%boundaryName == "none" )     THEN
                     this%nodes(this%edges(j)%nodeIds(k))%boundaryName = boundaryName
                  ELSE IF ( this%nodes(this%edges(j)%nodeIds(k))%boundaryName /= boundaryName )     THEN
                     this%nodes(this%edges(j)%nodeIds(k))%boundaryName = &
                          this%nodes(this%edges(j)%nodeIds(k))%boundaryName // boundaryName
                  END IF
               END DO
            ELSE
               this%edges(j)%edgeType = QMESH_INTERIOR
            END IF
         END DO

      END SUBROUTINE SetBoundaryTypes
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE FinishEdges( this, N ) 
         IMPLICIT NONE 
         TYPE(QuadMesh) :: this
         INTEGER        :: N
         INTEGER        :: j
!
!        ----------------------------------
!        Set loop indices for neighbor edge
!        ----------------------------------
!
         DO j = 1, this%numberOfEdges 
            IF( this%edges(j)%elementIDs(2) /= QMESH_NONE ) THEN
               IF ( this%edges(j)%elementSide(2) > 0 )     THEN
                  this%edges(j)%nStart = 1
                  this%edges(j)%nEnd   = N-1
                  this%edges(j)%nInc   = 1
               ELSE               
                  this%edges(j)%nStart = N-1
                  this%edges(j)%nEnd   = 1
                  this%edges(j)%nInc   = -1
               END IF
            END IF 
         END DO
      END SUBROUTINE FinishEdges
!
!     -------------
!     Print methods
!     -------------
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PrintMesh( this ) 
      IMPLICIT NONE 
      TYPE(QuadMesh) :: this
      INTEGER ::  k
      
      PRINT *, "Nodes..."
      DO k = 1, SIZE(this%nodes)
         CALL PrintNode( this%nodes(k), k )
      END DO
      PRINT *, "Elements..."
      DO k = 1, SIZE(this%elements) 
         CALL PrintElement( this%elements(k), k )
      END DO
      PRINT *, "Edges...s"
      DO k = 1, this%numberOfEdges 
         CALL PrintEdge( this%edges(k) )
      END DO
      
      END SUBROUTINE PrintMesh
      
!     ==========
      END MODULE QuadMeshClass
!     ==========
!
      INTEGER FUNCTION Hash1( idPair )
         INTEGER, DIMENSION(2) :: idPair
         Hash1 = MAXVAL(idPair)
      END FUNCTION Hash1
!
!////////////////////////////////////////////////////////////////////////
!
      INTEGER FUNCTION Hash2( idPair )
         INTEGER, DIMENSION(2) :: idPair
         Hash2 = MINVAL(idPair)
      END FUNCTION Hash2
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeLinearInterpolant( xStart, xEnd, bCurveOrder, nodes, values )
      USE SMConstants
      IMPLICIT NONE
      INTEGER       :: bCurveOrder
      REAL(KIND=RP) :: xStart(2), xEnd(2)
      REAL(KIND=RP) :: nodes(0:bCurveOrder), values(0:bCurveOrder,2)
      INTEGER       :: j
      
      DO j = 0, bCurveOrder 
         values(j,:) = xStart + (xEnd-xStart)*(nodes(j)+1.0_RP)*0.5_RP
      END DO
      END SUBROUTINE 
      