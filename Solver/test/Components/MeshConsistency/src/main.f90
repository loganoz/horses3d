!
!//////////////////////////////////////////////////////
!
!   @File:    main.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Dec 24 15:18:48 2017
!   @Last revision date: Tue Sep 25 10:03:25 2018
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: f80714f4c78c40aeae2f260a53c2dd6437a97038
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
   use BoundaryConditions,    only: BCs
   use FreeSlipWallBCClass,   only: FreeSlipWallBC_t
   use InterpolationMatrices, only: Initialize_InterpolationMatrices, Finalize_InterpolationMatrices
   use PhysicsStorage, only: SetReferenceLength
   implicit none
   TYPE(TestSuiteManager) :: testSuite
   INTEGER                :: numberOfFailures
   integer                :: i, k
   integer, parameter     :: Nmax = 40

   call MPI_Process % Init
   
   CALL testSuite % init()
   
   call SetReferenceLength(1.0_RP)
   call InitializeNodalStorage(GAUSS,Nmax)
   call Initialize_InterpolationMatrices(Nmax)
   CALL ConstructSharedBCModule
   allocate(BCs(12))
   do i = 1, 12
      allocate(FreeSlipWallBC_t    :: BCs(i) % bc) 
      select type(bc => BCs(i) % bc)
      type is (FreeSlipWallBC_t)
      bc % BCType = "freeslipwall"
      bc % isAdiabatic = .true.
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
   
   IF(numberOfFailures > 0)   error stop 99

end program main