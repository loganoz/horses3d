!
!////////////////////////////////////////////////////////////////////////
!
!      MeshFileTest.f90
!      Created: June 1, 2015 at 11:29 AM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE MeshFileTest  
         IMPLICIT NONE
         
         CALL readMeshFilewithName("SimpleBox.mesh")         
!         CALL readMeshFilewithName("BoxAroundCircle3D.mesh")
      END SUBROUTINE MeshFileTest
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE readMeshFilewithName(meshfileName)
         USE FTAssertions
         USE HexMeshClass 
         USE NodalStorageClass
         use ElementClass
         use ReadMeshFile
         IMPLICIT NONE  
         
         TYPE(HexMesh)                      :: mesh
         INTEGER                            :: N(3)
         INTEGER, ALLOCATABLE               :: Nvector(:)
         INTEGER                            :: nelem
         INTEGER                            :: fUnit
         INTEGER                            :: eID, l, NDOF, firstIdx
         CHARACTER(LEN=*)                   :: meshFileName
         LOGICAL                            :: success
         
         N = 6
         
         OPEN(newunit = fUnit, FILE = meshFileName )  
            READ(fUnit,*) l, nelem, l                    ! Here l is used as default reader since this variables are not important now
         CLOSE(fUnit)
         
         ALLOCATE (Nvector(nelem))
         Nvector = N(1)             ! No anisotropy
         call InitializeNodalStorage(GAUSS,N(1))
         
         call NodalStorage(N(1)) % Construct(GAUSS, N(1))
         call NodalStorage(N(2)) % Construct(GAUSS, N(2))
         call NodalStorage(N(3)) % Construct(GAUSS, N(3))
         CALL constructMeshFromFile(mesh,meshfileName,GAUSS, Nvector,Nvector,Nvector, .TRUE., 0, success)
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
      
!        Construct global storage
!        ------------------------
         call mesh % storage % construct(NDOF, 0)
      
!        Construct element storage
!        -------------------------
         firstIdx = 1
         DO eID = 1, SIZE(mesh % elements)
            associate (e => mesh % elements(eID))
            call mesh % elements(eID) % Storage % Construct(Nx = e % Nxyz(1), &
                                                            Ny = e % Nxyz(2), &
                                                            Nz = e % Nxyz(3), &
                                              computeGradients = .true., &
                                                 globalStorage = mesh % storage, &
                                                      firstIdx = firstIdx)
            firstIdx = firstIdx + e % Storage % NDOF
            end associate
         END DO
         
      END SUBROUTINE readMeshFilewithName
