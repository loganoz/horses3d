!
!//////////////////////////////////////////////////////
!
!   @File:    main.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Dec 24 15:18:48 2017
!   @Last revision date: Thu Jul 26 22:00:47 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: fb773e7c8706f4b4ef1f5bf9693a2b44f6c12dd2
!
!//////////////////////////////////////////////////////
!
program main
   use SMConstants
   use MeshTests
   use MeshConsistencySetup
   USE TestSuiteManagerClass
   use NodalStorageClass
   use SharedBCModule
   use MPI_Process_Info
   use BoundaryConditions, only: BCs
   use FreeSlipWallBCClass,   only: FreeSlipWallBC_t
   implicit none
   TYPE(TestSuiteManager) :: testSuite
   INTEGER                :: numberOfFailures
   integer                :: i
   integer, parameter     :: Nmax = 40

   call MPI_Process % Init
   
   CALL testSuite % init()
   
   call InitializeNodalStorage(Nmax)
   CALL ConstructSharedBCModule
   allocate(BCs(12))
   do i = 1, 12
      allocate(FreeSlipWallBC_t    :: BCs(i) % bc) 
      select type(bc => BCs(i) % bc)
      type is (FreeSlipWallBC_t)
      bc % BCType = "freeslipwall"
      bc % isAdiabatic = .true.
      bc % kWallType = 0.0_RP
      end select
   end do

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
