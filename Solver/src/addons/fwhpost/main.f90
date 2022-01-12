!This is file : main
! Author= oscarm
! Started at: 20.05.2021
! Last Modified: 20.05.2021
!
Program main
    use SMConstants
    use FTValueDictionaryClass
    use FileReaders               , only: ReadControlFile 
    use Headers
    use PhysicsStorage
    use DGSEMClass
    use StopwatchClass
    use pAdaptationClass          , only: GetMeshPolynomialOrders
    use NodalStorageClass
    use FluidData
    use InterpolationMatrices     , only: Initialize_InterpolationMatrices, Finalize_InterpolationMatrices
    use FWHPostProc
    use FWHPreSurface
    use SharedBCModule
    use MPI_Process_Info
#ifdef _HAS_MPI_
      use mpi
#endif
    Implicit None

    TYPE( FTValueDictionary)                :: controlVariables
    TYPE( DGSem )                           :: sem
    logical                                 :: success
    integer, allocatable                    :: Nx(:), Ny(:), Nz(:)
    integer                                 :: Nmax
    logical                                 :: surface

!   -----------------------------------------
!   Start measuring the total simulation time
!   -----------------------------------------
!
    call Stopwatch % CreateNewEvent("TotalTime")
    call Stopwatch % Start("TotalTime")

!   ---------------
!   Initializations
!   ---------------
!
    call MPI_Process % Init
!
!   ----------------------------------------------------------------------------------
!   The main is always compiled, so that __DATE__ and __TIME__ are updated accordingly
!   ----------------------------------------------------------------------------------
    call Main_Header("HORSES Ffowcs Williams Hawckings Analogy Post Proccesing",__DATE__,__TIME__)

    call controlVariables % initWithSize(16)
      CALL ConstructSharedBCModule

    call ReadControlFile( controlVariables )
    ! call CheckInputIntegrity(controlVariables, success)
    ! IF(.NOT. success)   ERROR STOP "Control file reading error"

!   ----------------
!   Set up the DGSEM
!   ----------------
!  
    ! Initialize manufactured solutions if necessary
    sem % ManufacturedSol = .FALSE.
      
    call ConstructPhysicsStorage( controlVariables, success )
    IF(.NOT. success)   ERROR STOP "Physics parameters input error"

    call GetMeshPolynomialOrders(controlVariables,Nx,Ny,Nz,Nmax)
    call InitializeNodalStorage (controlVariables ,Nmax)
    call Initialize_InterpolationMatrices(Nmax)

    call sem % construct (  controlVariables  = controlVariables,                                         &
        Nx_ = Nx,     Ny_ = Ny,     Nz_ = Nz,                                                 &
        success           = success)
    IF(.NOT. success)   ERROR STOP "Mesh reading error"

    surface = controlVariables % logicalValueForKey("fwhsurf")

    if (surface) then
        call extractSurface(sem % mesh, controlVariables)
    else
        call LoadAllFiles(controlVariables, sem)
    end if

!
!   ------------------------------------------
!   Finish measuring the total simulation time
!   ------------------------------------------
!
    call Stopwatch % Pause("TotalTime")
!
!   --------------------------
!   Show simulation statistics
!   --------------------------
!
   call DisplaySimulationStatistics(sem % numberOftimeSteps, sem % mesh)
!
!   ---------
!   Finish up
!   ---------
!
    call Stopwatch % destruct
    call sem % destruct()
    call Finalize_InterpolationMatrices
    call DestructGlobalNodalStorage()
   
End Program main

!
!//////////////////////////////////////////////////////////////////////// 
! 
      subroutine DisplaySimulationStatistics(iter,mesh)
         use SMConstants
         use HexMeshClass
         use StopwatchClass
         use Headers
         use MPI_Process_Info
#ifdef _HAS_MPI_
         use mpi
#endif
         implicit none
         integer,    intent(in)      :: iter
         type(HexMesh),   intent(in) :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: eID
         integer                    :: NDOF, localNDOF, ierr
         real(kind=RP)              :: Naverage, localNaverage
         real(kind=RP)              :: t_elaps, t_cpu
         
         if ( MPI_Process % isRoot ) write(STD_OUT,'(/)')
         call Section_Header("Simulation statistics")
         if ( MPI_Process % isRoot ) write(STD_OUT,'(/)')
!
!        Get mesh-related quantities
!        ---------------------------
         NDOF = 0
         Naverage = 0
   
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            NDOF = NDOF + (e % Nxyz(1) + 1)*(e % Nxyz(2) + 1)*(e % Nxyz(3) + 1)      
            Naverage = Naverage + e % Nxyz(1) + e % Nxyz(2) + e % Nxyz(3)
            end associate
         end do

         Naverage = Naverage / (3.0_RP * mesh % no_of_elements)
!
!        Perform a broadcast for the MPI solver
!        --------------------------------------
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
            localNDOF = NDOF
            localNaverage = Naverage * 3.0_RP * mesh % no_of_elements
            call mpi_allreduce(localNDOF, NDOF, 1, MPI_INT, MPI_SUM, &
                               MPI_COMM_WORLD, ierr)

            call mpi_allreduce(localNaverage, Naverage, 1, MPI_DOUBLE, MPI_SUM, &
                               MPI_COMM_WORLD, ierr)

            Naverage = Naverage / (3.0_RP * mesh % no_of_allElements)

         end if
#endif

         if ( .not. MPI_Process % isRoot ) return
!
!        Show preprocessing time
!        -----------------------
         t_elaps = Stopwatch % Elapsedtime("Preprocessing")
         t_cpu   = Stopwatch % CPUTime("Preprocessing")

         call Subsection_Header("Preprocessing")

         write(STD_OUT,'(30X,A,I0,A,F5.2,A,I0,A)')      "->   ", mesh % no_of_elements, &
                                                      " elements with polynomial order ",Naverage," (NDOF = ",NDOF,")."
         write(STD_OUT,'(30X,A,A30,ES10.3,A,ES10.3,A)') "->", "Preprocessing time: ",t_elaps," seconds (total CPU time: ",t_cpu,")."

!
!        Show simulation time
!        --------------------
         write(STD_OUT,'(/)')
         call Subsection_Header("TotalTime")
         if ( iter .le. 0 ) return

         t_elaps = Stopwatch % ElapsedTime("TotalTime")
         t_cpu   = Stopwatch % CPUTime("TotalTime")

         write(STD_OUT,'(30X,A,A30,ES10.3,A)') "->", "Simulation elapsed time: ",t_elaps," seconds."
         write(STD_OUT,'(30X,A,A30,ES10.3,A,ES10.3,A)') "->", "Simulation CPU time: ",t_cpu," seconds (ratio is ",t_cpu/t_elaps ,")."

      end subroutine DisplaySimulationStatistics
