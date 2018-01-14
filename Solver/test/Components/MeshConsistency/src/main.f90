!
!//////////////////////////////////////////////////////
!
!   @File:    main.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Dec 24 15:18:48 2017
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
program main
   use MeshTests
   use MeshConsistencySetup
   USE TestSuiteManagerClass
   use NodalStorageClass
   use SharedBCModule
   use MPI_Process_Info
   implicit none
   TYPE(TestSuiteManager) :: testSuite
   INTEGER                :: numberOfFailures
   integer                :: i
   integer, parameter     :: Nmax = 40

   call MPI_Process % Init
   
   CALL testSuite % init()
   
   call InitializeNodalStorage(Nmax)
   CALL ConstructSharedBCModule
   do i = 1, no_of_meshFiles
      call testSuite % addTestSubroutineWithName(OpenNextMesh,"Open next mesh: " // trim(meshFileNames(i)))
      call testSuite % addTestSubroutineWithName(CheckVolume,"Mesh volume: " // trim(meshFileNames(i)))
      call testSuite % addTestSubroutineWithName(CheckZoneSurfaces,"Zone surfaces: " // trim(meshFileNames(i)))
      call testSuite % addTestSubroutineWithName(CheckCoordinatesConsistency,&
                  "Mapping coordinates interpolation consistency: " // trim(meshFileNames(i)))
      call testSuite % addTestSubroutineWithName(CheckNormalConsistency,&
                  "Consistency in the normal vectors: " // trim(meshFileNames(i)))
      call testSuite % addTestSubroutineWithName(CheckScalConsistency,&
                  "Consistency in the surface Jacobians: " // trim(meshFileNames(i)))

   end do

   CALL testSuite % performTests(numberOfFailures)
   CALL testSuite % finalize()
   
   IF(numberOfFailures > 0)   STOP 99

end program main
