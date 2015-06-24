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
         USE PlotterClass
         USE NodalStorageClass
         IMPLICIT NONE  
         
         TYPE(NodalStorage)                 :: spA
         TYPE(HexMesh)                      :: mesh
         TYPE(Plotter)                      :: meshPlotter
         CLASS(PlotterDataSource), POINTER  :: dataSource
         INTEGER                            :: N = 6
         INTEGER                            :: id, l
         CHARACTER(LEN=*)                   :: meshFileName
         CHARACTER(LEN=128)                 :: plotFileName
         LOGICAL                            :: success
         
         CALL ConstructNodalStorage(spA, N)
         CALL mesh % constructFromFile(meshfileName,spA, success)
         CALL FTAssert(test = success,msg = "Mesh file read properly")
         IF(.NOT. success) RETURN 
         
         DO id = 1, SIZE(mesh % elements)
            CALL allocateElementStorage(this = mesh % elements(id),&
                                        N = N, nEqn = 3,nGradEqn = 0,flowIsNavierStokes = .FALSE.) 
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
         CALL meshPlotter % Construct(fUnit = 11,dataSource = dataSource)
         CALL meshPlotter % ExportToTecplot(elements = mesh % elements)
         CALL meshPlotter % Destruct()
         DEALLOCATE(dataSource)
         
      END SUBROUTINE readMeshFilewithName
