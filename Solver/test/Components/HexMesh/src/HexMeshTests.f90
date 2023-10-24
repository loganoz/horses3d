!
!////////////////////////////////////////////////////////////////////////
!
!      HexMeshTests.f90
!      Created: May 27, 2015 at 2:44 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
! 
      SUBROUTINE testTwoBoxesMeshConstruction 
         USE FTAssertions
         USE SMConstants
         use FaceClass
         USE HexMeshClass 
         use MeshTypes
         USE SharedBCModule
         use ReadMeshFile
         use TransfiniteMapClass
         use NodalStorageClass
         use ElementClass 
         use BoundaryConditions
         use FreeSlipWallBCClass
         use InterpolationMatrices, only: Initialize_InterpolationMatrices, Finalize_InterpolationMatrices
         IMPLICIT NONE
         
         TYPE(HexMesh), target           :: mesh
         TYPE(TransfiniteHexMap)         :: hexTransform
         TYPE(Face)                      :: testFace
         REAL(KIND=RP)                   :: nodes(3,12)
         REAL(KIND=RP)                   :: corners(3,8)
         REAL(KIND=RP)                   :: erMax
         REAL(KIND=RP)                   :: x(3), p(3)
         REAL(KIND=RP)                   :: u(3)
         INTEGER                         :: nodeIDs(8,2)
         INTEGER                         :: nonne(6)
         INTEGER                         :: i,j,k, N(3), l, id
         INTEGER                 :: nElements, nNodes, nFaces, iFaceID
         INTEGER                            :: eID, NDOF, firstIdx
         INTEGER                 :: numberOfBoundaryFaces, NumberofInteriorFaces
         CHARACTER(LEN=1)        :: space = " "
         CHARACTER(LEN=16)       :: meshfileName = "TwoElements.mesh"
         CHARACTER(LEN=128)      :: msg = "Node id xxx"
         LOGICAL                 :: success
         
         CHARACTER(len=8)                   :: boundaryNames(6)
         INTEGER, ALLOCATABLE               :: Nvector(:)
         INTEGER                            :: nelem
         INTEGER                            :: fUnit
         
         N           = 6
!
!        --------------------------------------
!        Create a simple mesh with two elements
!        --------------------------------------
!
         nodes(:,1)  = [0.0_RP, 0.0_RP, 0.0_RP]
         nodes(:,2)  = [1.0_RP, 0.0_RP, 0.0_RP]
         nodes(:,3)  = [1.0_RP, 2.0_RP, 0.0_RP]
         nodes(:,4)  = [0.0_RP, 2.0_RP, 0.0_RP]
         nodes(:,5)  = [0.0_RP, 0.0_RP, 3.0_RP]
         nodes(:,6)  = [1.0_RP, 0.0_RP, 3.0_RP]
         nodes(:,7)  = [1.0_RP, 2.0_RP, 3.0_RP]
         nodes(:,8)  = [0.0_RP, 2.0_RP, 3.0_RP]
         nodes(:,9)  = [2.0_RP, 0.0_RP, 0.0_RP]
         nodes(:,10) = [2.0_RP, 2.0_RP, 0.0_RP]
         nodes(:,11) = [2.0_RP, 0.0_RP, 3.0_RP]
         nodes(:,12) = [2.0_RP, 2.0_RP, 3.0_RP]
         
         nodeIDs(:,1) = [1,2,3,4,5,6,7,8]
         nodeIDs(:,2) = [2,9,10,3,6,11,12,7]
         nonne        = 0
         
         OPEN(UNIT = 11, FILE = meshfileName)
!
!        -------------------------
!        Write out the mesh header
!        -------------------------
!
         WRITE(11,*) 12, 2, N(1)
!
!        -------------------
!        Write out the nodes
!        -------------------
!
         DO i = 1, 12
            WRITE(11,*) nodes(:,i) 
         END DO  
!
!        ----------------------
!        Write out the elements
!
!        Elements have the format:
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
         WRITE(11,*) nodeIDs(:,1)
         WRITE(11,*) nonne
         WRITE(11,*) "front", space,  "back", space, "bottom", space, "---", space, "top", space, "left"
         
         WRITE(11,*) nodeIDs(:,2)
         WRITE(11,*) nonne
         WRITE(11,*) "front", space, "back", space, "bottom", space, "right", space, "top", space, "---"
         CLOSE(11)
!
!        -----------------
!        Generate the mesh
!        -----------------
!
         ! First, we have to generate the bcTypeDictionary as it would be read from the control file
         call constructSharedBCModule()
         allocate(BCs(6))
         do id = 1, 6
            allocate(FreeSlipWallBC_t    :: BCs(id) % bc) 
            select type(bc => BCs(id) % bc)
            type is (FreeSlipWallBC_t)
            bc % BCType = "freeslipwall"
            bc % isAdiabatic = .true.
            end select
         end do
         
         OPEN(newunit = fUnit, FILE = meshFileName )  
            READ(fUnit,*) l, nelem, l                    ! Here l is used as default reader since this variables are not important now
         CLOSE(fUnit)
         
         ALLOCATE (Nvector(nelem))
         Nvector = N(1)
         call InitializeNodalStorage(GAUSS, N(1))
         call Initialize_InterpolationMatrices(N(1))
         
         call NodalStorage(N(1)) % Construct(GAUSS, N(1))
         call NodalStorage(N(2)) % Construct(GAUSS, N(2))
         call NodalStorage(N(3)) % Construct(GAUSS, N(3))
         CALL constructMeshFromFile(mesh, meshfileName,GAUSS,Nvector,Nvector,Nvector, .TRUE.,  0, .false., success)
         
         CALL FTAssert(test = success,msg = "Mesh file properly constructed")
         IF(.NOT. success) return

!
!        ****************
!        Allocate storage: since there are no control variables, I have copied the AllocateStorage function
!        ****************
!
         NDOF = 0
         do eID = 1, size(mesh % elements)
            associate (e => mesh % elements(eID))
            NDOF = NDOF + (e % Nxyz(1) + 1)*(e % Nxyz(2) + 1)*(e % Nxyz(3) + 1) 
            end associate
         end do
!!      
!!        Construct global storage
!!        ------------------------
!         call mesh % storage % construct (NDOF, Nvector, Nvector, Nvector, .TRUE., .FALSE. )
!         DO eID = 1, nelem
!            mesh % elements(eID) % storage => mesh % storage % elements(eID)
!         END DO
!
!         call mesh % SetStorageToEqn(1)
!
!        ---------------------------
!        Check integrity of the mesh
!        ---------------------------
!
         nNodes = SIZE( mesh % nodes )
         CALL FTAssertEqual(expectedValue = 12 ,actualValue = nNodes, msg = "Number of nodes in mesh")
         nElements = SIZE( mesh % elements)
         CALL FTAssertEqual(expectedValue = 2 ,actualValue = nElements, msg = "Number of elements in mesh")
         
         CALL FTAssertEqual(expectedValue = NodalStorage(N(1)) % N, &
                                                  actualValue = mesh % elements(1) % Nxyz(1), msg = "Polynomial order of element 1")
         CALL FTAssertEqual(expectedValue = NodalStorage(N(1)) % N, &
                                                  actualValue = mesh % elements(2) % Nxyz(1), msg = "Polynomial order of element 2")
         
         DO k = 1, 2
            DO i = 1, 8
               WRITE(msg,*) "Node id ",i, ", element", k
               CALL FTAssertEqual(expectedValue = nodeIDs(i,k),actualValue = mesh % elements(k) % nodeIDs(i), msg = TRIM(msg)) 
            END DO  
         END DO
         
         nFaces = SIZE(mesh % faces)
         CALL FTAssertEqual(expectedValue = 11 ,actualValue = nFaces, msg = "Number of faces in mesh")
!
!        ----------------------------
!        Check the element geometries
!        ----------------------------
!
         erMax = 0.0_RP
         DO id = 1, 2
            N = mesh % elements(id) % Nxyz
            DO l = 1, 8
               corners(:,l) = nodes(:,nodeIDs(l,id))
            END DO   
            CALL hexTransform % constructWithCorners(corners)
            
            DO k = 0, mesh % elements(id) % Nxyz(3)
               u(3) = NodalStorage(N(3)) % x(k)
               DO j = 0, mesh % elements(id) % Nxyz(2)
                  u(2) = NodalStorage(N(2)) % x(j)
                  DO i = 0, mesh % elements(id) % Nxyz(1)
                     u(1) = NodalStorage(N(1)) % x(i)
                     p = hexTransform % transfiniteMapAt(u)
                     x = mesh % elements(id) % geom % x(:,i,j,k)
                     erMax = MAX(erMax,MAXVAL(ABS(x-p)))
                  END DO
               END DO
            END DO  
            
            WRITE(msg,*) "element collocation locations in element ",id
            CALL FTAssertEqual(expectedValue = 0.0_RP, actualValue = erMax, tol = 1.d-10, msg = msg)
            CALL hexTransform % destruct()
            
         END DO
!
!        -----------
!        Check faces
!        -----------
!
         CALL FTAssertEqual(expectedValue = 11 ,actualValue = mesh % numberOfFaces, msg = "Number of faces in mesh")
         
         NumberofInteriorFaces = 0
         numberOfBoundaryFaces = 0
         DO j = 1, SIZE(mesh % faces)
            IF( mesh % faces(j) % elementIDs(2) == HMESH_NONE)     THEN
               numberOfBoundaryFaces = numberOfBoundaryFaces + 1
            ELSE
               NumberofInteriorFaces = NumberofInteriorFaces + 1
               iFaceID = j
            END IF  
         END DO  
         CALL FTAssertEqual(expectedValue = 10 ,actualValue = numberOfBoundaryFaces, msg = "Number of boundary faces in mesh")
         CALL FTAssertEqual(expectedValue = 1  ,actualValue = NumberofInteriorFaces, msg = "Number of interior faces in mesh")
         
          
         testFace = mesh % faces(iFaceID)
        
         CALL FTAssertEqual(expectedValue = 0, actualValue = testFace % rotation,msg = "rotation of slave to master face")
         CALL FTAssertEqual(expectedValue = 1, actualValue = testFace % elementIDs(1),msg = "Face element master")
         CALL FTAssertEqual(expectedValue = 2, actualValue = testFace % elementIDs(2),msg = "Face element slave")
         CALL FTAssertEqual(expectedValue = 4, actualValue = testFace % elementSide(1),msg = "Face element master side")
         CALL FTAssertEqual(expectedValue = 6, actualValue = testFace % elementSide(2),msg = "Face element slave side")
!
!        ---------
!        Finish up
!        ---------
!
         OPEN(UNIT=11, FILE = meshfileName)
         CLOSE(11, STATUS = "DELETE")
         
!
!        --------------------------------------------------------------
!        Call to DestructGlobalNodalStorage created an error with ifort
!        Copy the subroutine here to solved the issue         
!        --------------------------------------------------------------
!

         !call DestructGlobalNodalStorage()

         if ( allocated(NodalStorage_Gauss) ) then
            do k=lbound(NodalStorage_Gauss,1), ubound(NodalStorage_Gauss,1)
               IF (.NOT. NodalStorage_Gauss(k) % Constructed) cycle
               call NodalStorage_Gauss(k) % destruct()
            end do
            deallocate (NodalStorage_Gauss)
         end if
   
         if ( allocated(NodalStorage_GaussLobatto) ) then
            do k=lbound(NodalStorage_GaussLobatto,1), ubound(NodalStorage_GaussLobatto,1)
               IF (.NOT. NodalStorage_GaussLobatto(k) % Constructed) cycle
               call NodalStorage_GaussLobatto(k) % destruct()
            end do      
            deallocate (NodalStorage_GaussLobatto)
         end if
   
         nullify (NodalStorage)

         call Finalize_InterpolationMatrices
         
      END SUBROUTINE testTwoBoxesMeshConstruction
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE testTwoElementCylindersMesh  
         USE FTAssertions
         USE SMConstants
         USE HexMeshClass 
         use ReadMeshFile
         use FaceClass
         use NodalStorageClass
         use InterpolationMatrices, only: Initialize_InterpolationMatrices, Finalize_InterpolationMatrices
         use ElementClass
         use MeshTypes
         IMPLICIT NONE
         
         EXTERNAL                :: cylindricalGeometry
         TYPE(HexMesh), target   :: mesh
         TYPE(Face)              :: testFace
         INTEGER                 :: j, k, N(3), id
         INTEGER                 :: iFaceID
         INTEGER                 :: numberOfBoundaryFaces, NumberofInteriorFaces
         INTEGER                 :: eID, l, NDOF, firstIdx
         CHARACTER(LEN=27)       :: meshfileName = "TwoClyindricalElements.mesh"
         LOGICAL                 :: success
         
         INTEGER, ALLOCATABLE               :: Nvector(:)
         INTEGER                            :: nelem
         INTEGER                            :: fUnit
!         
         N = 6
!
!        ------------------
!        Create a mesh file
!        ------------------
!
         CALL generateClyindricalMeshFile(meshfileName)
!
!        -------------------------------
!        Generate the mesh from the file
!        -------------------------------
!
         OPEN(newunit = fUnit, FILE = meshFileName )  
            READ(fUnit,*) l, nelem, l                    ! Here l is used as default reader since this variables are not important now
         CLOSE(fUnit)
         
         ALLOCATE (Nvector(nelem))
         Nvector = N(1)
         call InitializeNodalStorage(GAUSS, N(1))
         call Initialize_InterpolationMatrices(N(1))
         
         call NodalStorage(N(1)) % Construct(GAUSS, N(1))
         call NodalStorage(N(2)) % Construct(GAUSS, N(2))
         call NodalStorage(N(3)) % Construct(GAUSS, N(3))
         CALL constructMeshFromFile(mesh,meshfileName,GAUSS,Nvector,Nvector,Nvector, .TRUE., 0, .false., success)
         
         CALL FTAssert(test = success,msg = "Mesh file read properly")
         IF(.NOT. success) RETURN 
!
!        ****************
!        Allocate storage: since there are no control variables, I have copied the AllocateStorage function
!        ****************
!
         NDOF = 0
         do eID = 1, size(mesh % elements)
            associate (e => mesh % elements(eID))
            NDOF = NDOF + (e % Nxyz(1) + 1)*(e % Nxyz(2) + 1)*(e % Nxyz(3) + 1) 
            end associate
         end do
         
         call mesh % storage % construct (NDOF, Nvector, Nvector, Nvector, .TRUE., .FALSE. )
         DO eID = 1, nelem
            mesh % elements(eID) % storage => mesh % storage % elements(eID)
         END DO

         CALL FTAssertEqual(expectedValue = 2     ,actualValue = SIZE(mesh % elements),msg = "Number of elements in mesh.")
         CALL FTAssertEqual(expectedValue = N(1)+1,actualValue = SIZE(mesh % elements(1) % storage % Q,2) ,msg = "Number of solution points")
!
!        -----------
!        Check faces
!        -----------
!
         CALL FTAssertEqual(expectedValue = 11 ,actualValue = mesh % numberOfFaces, msg = "Number of faces in mesh")
         
         NumberofInteriorFaces = 0
         numberOfBoundaryFaces = 0
         DO j = 1, SIZE(mesh % faces)
            IF( mesh % faces(j) % elementIDs(2) == HMESH_NONE)     THEN
               numberOfBoundaryFaces = numberOfBoundaryFaces + 1
            ELSE
               NumberofInteriorFaces = NumberofInteriorFaces + 1
               iFaceID = j
            END IF  
         END DO  
         CALL FTAssertEqual(expectedValue = 10 ,actualValue = numberOfBoundaryFaces, msg = "Number of boundary faces in mesh")
         CALL FTAssertEqual(expectedValue = 1  ,actualValue = NumberofInteriorFaces, msg = "Number of interior faces in mesh")
         
          
         testFace = mesh % faces(iFaceID)
        
         CALL FTAssertEqual(expectedValue = 0, actualValue = testFace % rotation     ,msg  = "rotation of slave to master face")
         CALL FTAssertEqual(expectedValue = 1, actualValue = testFace % elementIDs(1),msg  = "Face element master")
         CALL FTAssertEqual(expectedValue = 2, actualValue = testFace % elementIDs(2),msg  = "Face element slave")
         CALL FTAssertEqual(expectedValue = 5, actualValue = testFace % elementSide(1),msg = "Face element master side")
         CALL FTAssertEqual(expectedValue = 3, actualValue = testFace % elementSide(2),msg = "Face element slave side")
         
!
!        --------------------------------------------------------------
!        Call to DestructGlobalNodalStorage created an error with ifort
!        Copy the subroutine here to solved the issue         
!        --------------------------------------------------------------
!

         !call DestructGlobalNodalStorage()

         if ( allocated(NodalStorage_Gauss) ) then
            do k=lbound(NodalStorage_Gauss,1), ubound(NodalStorage_Gauss,1)
               IF (.NOT. NodalStorage_Gauss(k) % Constructed) cycle
               call NodalStorage_Gauss(k) % destruct()
            end do
            deallocate (NodalStorage_Gauss)
         end if
   
         if ( allocated(NodalStorage_GaussLobatto) ) then
            do k=lbound(NodalStorage_GaussLobatto,1), ubound(NodalStorage_GaussLobatto,1)
               IF (.NOT. NodalStorage_GaussLobatto(k) % Constructed) cycle
               call NodalStorage_GaussLobatto(k) % destruct()
            end do      
            deallocate (NodalStorage_GaussLobatto)
         end if
   
         nullify (NodalStorage)

         call Finalize_InterpolationMatrices
      END SUBROUTINE testTwoElementCylindersMesh
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE generateClyindricalMeshFile(meshfileName)
         USE SMConstants
         USE FacePatchClass
         IMPLICIT NONE
         
         CHARACTER(LEN=*)        :: meshfileName
         EXTERNAL                :: cylindricalGeometry
         REAL(KIND=RP)           :: nodes(3,12)
         INTEGER                 :: nodeIDs(8,2), curvedFaces(6)
         INTEGER                 :: i,j,k, N(3)
         CHARACTER(LEN=1)        :: space = " "
!
!        -----
!        Faces
!        -----
!
         INTEGER         , PARAMETER        :: LEFT_FACE = 6, RIGHT_FACE = 4, TOP_FACE = 5, BOTTOM_FACE = 3
         INTEGER                            :: nUknots, nVKnots
         REAL(KIND=RP)   , ALLOCATABLE      :: uKnots(:)
         REAL(KIND=RP)   , ALLOCATABLE      :: vKnots(:)
         REAL(KIND=RP)   , ALLOCATABLE      :: points(:,:,:)
         TYPE(FacePatch) , DIMENSION(6)     :: faceData
!         
         N       = 6
         nUknots = N(1) + 1
         nVKnots = N(2) + 1
         ALLOCATE(uKnots(nUKnots))
         ALLOCATE(vKnots(nVKnots))
         ALLOCATE(points(3,nUKnots,nVKnots))
         
         DO i = 1, nUKnots
            uKnots(i) = -COS((i-1)*PI/(nUKnots-1)) 
         END DO  
         
         DO j = 1, nVKnots
            vKnots(j) = -COS((j-1)*PI/(nVKnots-1)) 
         END DO  
!
!        -------------------------------------
!        Create four cylindrically curved faces
!        -------------------------------------
!
         DO j = 1, nVKnots
            DO i = 1, nUKnots
               CALL cylindricalGeometry([-1.0_RP, uKnots(i), vKnots(j)], points(:,i,j))
            END DO   
         END DO
         CALL ConstructFacePatch( faceData(LEFT_FACE), uKnots, vKnots, points )
         
         DO j = 1, nVKnots
            DO i = 1, nUKnots
               CALL cylindricalGeometry([1.0_RP, uKnots(i), vKnots(j)], points(:,i,j))
            END DO   
         END DO
         CALL ConstructFacePatch( faceData(RIGHT_FACE), uKnots, vKnots, points )
         
         DO j = 1, nVKnots
            DO i = 1, nUKnots
               CALL cylindricalGeometry([uKnots(i), vKnots(j), -1.0_RP], points(:,i,j))
            END DO   
         END DO
         CALL ConstructFacePatch( faceData(BOTTOM_FACE), uKnots, vKnots, points )
         
         DO j = 1, nVKnots
            DO i = 1, nUKnots
               CALL cylindricalGeometry([uKnots(i), vKnots(j), 1.0_RP ], points(:,i,j))
            END DO   
         END DO
         CALL ConstructFacePatch( faceData(TOP_FACE), uKnots, vKnots, points )
!
!        ----------------------
!        Write out the elements
!
!        Elements have the format:
!        node1ID node2ID node3ID node4ID ... node8ID
!        b1 b2 b3 b4 ... b8
!           (=0 for straight side, 1 for curved)
!        If curved boundaries, then for each:
!           for j = 1 to nFacePoints
!               for i = 1, nFacePoints
!                    x_(i,j)  y_(i,j) z_(i,jj)
!               next i
!           next j
!        bname1 bname2 bname3 bname4 ... bname8
!        -----------------------------------------
!
         nodes(:,1)  = faceData(LEFT_FACE)  % points(:,1,1)
         nodes(:,2)  = faceData(RIGHT_FACE) % points(:,1,1)
         nodes(:,3)  = faceData(RIGHT_FACE) % points(:,nUknots,1)
         nodes(:,4)  = faceData(LEFT_FACE)  % points(:,nUknots,1)
         nodes(:,5)  = faceData(LEFT_FACE)  % points(:,1,nUknots)
         nodes(:,6)  = faceData(RIGHT_FACE) % points(:,1,nUknots)
         nodes(:,7)  = faceData(RIGHT_FACE) % points(:,nUknots,nUknots)
         nodes(:,8)  = faceData(LEFT_FACE)  % points(:,nUknots,nUknots)
         nodes(:,9)  = nodes(:,5)
         nodes(:,10) = nodes(:,6)
         nodes(:,11) = nodes(:,7)
         nodes(:,12) = nodes(:,8)
         nodes(3,9:12) = nodes(3,9:12) + 3.0_RP
         
         nodeIDs(:,1) = [1,2,3,4,5,6,7,8]
         nodeIDs(:,2) = [5,6,7,8,9,10,11,12]
         curvedFaces  = [0,0,1,1,1,1]
         
         OPEN(UNIT = 11, FILE = meshfileName)
!
!        ------------------
!        Write out the mesh
!        ------------------
!
         WRITE(11,*) 12, 2, N(1)
!
!        -------------------
!        Write out the nodes
!        -------------------
!
         DO i = 1, 12
            WRITE(11,*) nodes(:,i) 
         END DO  
         
         WRITE(11,*) nodeIDs(:,1)
         WRITE(11,*) curvedFaces
         DO k = 1, 6
            IF ( curvedFaces(k) == 1 )     THEN
               DO j = 1, nVKnots
                  DO i = 1, nVKnots
                     WRITE(11,*) faceData(k) % points(:,i,j) 
                  END DO   
               END DO  
            END IF
         END DO  
         WRITE(11,*) "front", space,  "back", space, "bottom", space, "right", space, "---", space, "left"
         
         WRITE(11,*) nodeIDs(:,2)
         WRITE(11,*) curvedFaces
         
         DO k = 1, 6
            IF ( curvedFaces(k) == 1 )     THEN
               faceData(k) % points(3,:,:) = faceData(k) % points(3,:,:) + 3.0_RP !shift defined in CylindricalGeometry
            END IF  
         END DO  
         DO k = 1, 6
            IF ( curvedFaces(k) == 1 )     THEN
               DO j = 1, nVKnots
                  DO i = 1, nVKnots
                     WRITE(11,*) faceData(k) % points(:,i,j) 
                  END DO   
               END DO  
            END IF
         END DO  
         
         WRITE(11,*) "front", space, "back", space, "---", space, "right", space, "top", space, "left"
         CLOSE(11)
         
      END SUBROUTINE generateClyindricalMeshFile
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE cylindricalGeometry(u,x)
         USE SMConstants
         IMPLICIT NONE  
         REAL(KIND=RP) :: u(3) ! IN [-1,1]^3
         REAL(KIND=RP) :: x(3) ! physical space locations
         REAL(KIND=RP) :: r0 = 1.0_RP, rMax = 2.0_RP, theta0 = 0.0_RP, &
                          thetaMax = PI/2, z0 = 0.0_RP, zMax = 3.0_RP
         REAL(KIND=RP) :: r, theta, z
         
         r     = r0 + 0.5_RP*(u(1)+1.0_RP)*(rMax - r0)
         theta = theta0 + 0.5_RP*(u(2)+1.0_RP)*(thetaMax - theta0)
         z     = z0 + 0.5_RP*(u(3)+1.0_RP)*(zMax - z0)
         
         x(1) = r*COS(theta)
         x(2) = r*SIN(theta)
         x(3) = z
         
      END SUBROUTINE cylindricalGeometry