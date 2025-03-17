Module readSpecM
    use HexMeshClass
    use ElementClass
    use SMConstants
    use ElementConnectivityDefinitions
    use LocalRefinement
    use Utilities, only: UnusedUnit
    use PhysicsStorage
    use NodeClass
    Implicit None

    contains

    Subroutine ConstructSimpleMesh_FromSpecFile_(self, fileName, locR, AllNx, AllNy, AllNz)

!        ---------------
!        Input variables
!        ---------------
!
         type(HexMesh)                    :: self
         CHARACTER(LEN=*)                 :: fileName
         type(LocalRef_t), optional, intent(in)  :: locR
         integer         , optional, intent(in)  :: AllNx(:), AllNy(:), AllNz(:)     !<  Polynomial orders for all the elements
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
         ! CHARACTER(LEN=BC_STRING_LENGTH), pointer :: zoneNames(:)
         real(kind=RP)                   :: corners(NDIM,NODES_PER_ELEMENT)
         integer                         :: Nx, Ny, Nz     !<  Polynomial orders for each element
         integer, dimension(8)      :: falseNodeID          ! dummy variable needed to construct element, the actual nodes are not computed
!
!        ------------------
!        For curved patches
!        ------------------
!
         ! real(kind=RP)  , DIMENSION(:)    , ALLOCATABLE :: uNodes, vNodes
         real(kind=RP)  , DIMENSION(:,:,:), ALLOCATABLE :: values
!
!        ----------------
!        For flat patches
!        ----------------
!
         ! real(kind=RP)  , DIMENSION(2)     :: uNodesFlat = [-1.0_RP,1.0_RP]
         ! real(kind=RP)  , DIMENSION(2)     :: vNodesFlat = [-1.0_RP,1.0_RP]
         ! real(kind=RP)  , DIMENSION(3,2,2) :: valuesFlat


         numberOfBoundaryFaces = 0
!
!        -----------------------
!        Read header information
!        -----------------------
!
         fUnit = UnusedUnit()
         OPEN( UNIT = fUnit, FILE = fileName, iostat = fileStat )
         IF ( fileStat /= 0 )     THEN
            PRINT *, "Error opening file: ", fileName
            RETURN 
         END IF
         
         READ(fUnit,*) numberOfNodes, numberOfElements, bFaceOrder

         self % no_of_elements = numberOfElements
         self % no_of_allElements = numberOfElements

         numBFacePoints = bFaceOrder + 1
         ! allocate(uNodes(numBFacePoints))
         ! allocate(vNodes(numBFacePoints))
         allocate(values(3,numBFacePoints,numBFacePoints))
         
         ! DO i = 1, numBFacePoints
         !    uNodes(i) = -COS((i-1.0_RP)*PI/(numBFacePoints-1.0_RP)) 
         !    vNodes(i) = -COS((i-1.0_RP)*PI/(numBFacePoints-1.0_RP)) 
         ! END DO
!
!        ---------------
!        Allocate memory
!        ---------------
!
         allocate( self % elements(numberOfelements) )
         allocate( self % nodes(numberOfNodes) )
!
!        ----------------------------------
!        Read nodes: Nodes have the format:
!        x y z
!        ----------------------------------
!
         DO j = 1, numberOfNodes 
            READ( fUnit, * ) x
            ! x = x/Lref
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
               
            ELSE

            ! get corners for simple mesh to do the local refinement
               DO k = 1, NODES_PER_ELEMENT
                  corners(:,k) = self % nodes(nodeIDs(k)) % x
               END DO
!
!              --------------------------------------------------------------
!              Otherwise, we have to look at each of the faces of the element 
!              --------------------------------------------------------------
!
               DO k = 1, FACES_PER_ELEMENT 
                 IF ( faceFlags(k) == 0 )     THEN
!
               ! nothing to do here
                  ELSE
!
!                    -------------
!                    Curved faces 
!                    -------------
!
                    ! read values and do nothing with them, useful to go to next important value of file
                     DO j = 1, numBFacePoints
                        DO i = 1, numBFacePoints
                           READ( fUnit, * ) values(:,i,j)
                        END DO  
                     END DO
   
                  END IF
               END DO
               
            END IF 
!
!           -------------------------
!           Now construct the element
!           -------------------------
!
            ! set dummy values
            falseNodeID = 0

            if( present(locR) ) then 
               call locR % getOrderOfPosition(corners, Nx, Ny, Nz)
               ! call self % elements(l) % Construct (Nx, Ny, Nz, falseNodeID , l, l) 
               call self % elements(l) % Construct (Nx, Ny, Nz, nodeIDs, l, l) 
            else
               call self% elements(l)% Construct(AllNx(l), AllNy(l), AllNz(l), nodeIDs, l, l) 
            end if

            READ( fUnit, * ) names
        END DO

         CLOSE( fUnit )

    End Subroutine ConstructSimpleMesh_FromSpecFile_

End Module readSpecM
