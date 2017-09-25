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
         USE DGSEMPlotterClass
         USE NodalStorageClass
         IMPLICIT NONE  
         
         TYPE(NodalStorage), ALLOCATABLE    :: spA(:,:,:)
         TYPE(HexMesh)                      :: mesh
         TYPE(DGSEMPlotter)                 :: meshPlotter
         CLASS(PlotterDataSource), POINTER  :: dataSource
         INTEGER                            :: N(3)
         INTEGER, ALLOCATABLE               :: Nvector(:)
         INTEGER                            :: nelem
         INTEGER                            :: fUnit
         INTEGER                            :: id, l
         CHARACTER(LEN=*)                   :: meshFileName
         CHARACTER(LEN=128)                 :: plotFileName
         LOGICAL                            :: success
         
         N = 6
         
         ALLOCATE(spA(0:N(1),0:N(2),0:N(3)))
         OPEN(newunit = fUnit, FILE = meshFileName )  
            READ(fUnit,*) l, nelem, l                    ! Here l is used as default reader since this variables are not important now
         CLOSE(fUnit)
         
         ALLOCATE (Nvector(nelem))
         Nvector = N(1)             ! No anisotropy
         
         CALL ConstructNodalStorage(spA(N(1),N(2),N(3)), N(1),N(2),N(3))
         CALL mesh % constructFromFile(meshfileName,spA, Nvector,Nvector,Nvector, success)
         CALL FTAssert(test = success,msg = "Mesh file read properly")
         IF(.NOT. success) RETURN 
         
         DO id = 1, SIZE(mesh % elements)
            CALL allocateElementStorage(self = mesh % elements(id),&
                                        Nx = N(1), Ny = N(2), Nz = N(3), nEqn = 5,nGradEqn = 0,flowIsNavierStokes = .FALSE.) 
         END DO
         
!
!        ---------------------
!        Plot the mesh read in
!        ---------------------
!
         l = INDEX(STRING = meshFileName, SUBSTRING = ".")
         plotFileName = meshFileName(1:l) // "tec"
         ALLOCATE(dataSource)
         OPEN(UNIT = 11, FILE = plotFileName)
         CALL meshPlotter % Construct  (fUnit      = 11,          &
                                        dataSource = dataSource,      &
                                        newN       = N(1), &
                                        spA        = spA)
         
         
         CALL meshPlotter % ExportToTecplot(elements = mesh % elements, spA = spA)
         CALL meshPlotter % Destruct()
         DEALLOCATE(dataSource)
         
      END SUBROUTINE readMeshFilewithName
