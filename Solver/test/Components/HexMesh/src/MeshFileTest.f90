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
         use InterpolationMatrices, only: Initialize_InterpolationMatrices, Finalize_InterpolationMatrices
         IMPLICIT NONE  
         
         TYPE(HexMesh), target              :: mesh
         INTEGER                            :: N(3)
         INTEGER, ALLOCATABLE               :: Nvector(:)
         INTEGER                            :: nelem
         INTEGER                            :: fUnit
         INTEGER                            :: eID, l, NDOF, k
         CHARACTER(LEN=*)                   :: meshFileName
         LOGICAL                            :: success
         
         N = 6
         
         OPEN(newunit = fUnit, FILE = meshFileName )  
            READ(fUnit,*) l, nelem, l                    ! Here l is used as default reader since this variables are not important now
         CLOSE(fUnit)
         
         ALLOCATE (Nvector(nelem))
         Nvector = N(1)             ! No anisotropy
         call InitializeNodalStorage(GAUSS,N(1))
         call Initialize_InterpolationMatrices(N(1))
         
         call NodalStorage(N(1)) % Construct(GAUSS, N(1))
         call NodalStorage(N(2)) % Construct(GAUSS, N(2))
         call NodalStorage(N(3)) % Construct(GAUSS, N(3))
         CALL constructMeshFromFile(mesh,meshfileName,GAUSS, Nvector,Nvector,Nvector, .TRUE., 0, .false., success)
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
      END SUBROUTINE readMeshFilewithName