program horses2plt
   use SMConstants
   use getTask
   use Mesh2PltModule
   use Solution2PltModule
   use Stats2PltModule
   use SolutionFile
   use SharedSpectralBasis
   implicit none
   integer                                 :: jobType
   character(len=LINE_LENGTH)              :: meshName
   integer                                 :: no_of_solutions
   character(len=LINE_LENGTH), allocatable :: solutionNames(:)
   integer, allocatable                    :: solutionTypes(:)
   logical                                 :: fixedOrder
   integer                                 :: Nout(3)
   integer                                 :: basis
   integer                                 :: iSol
!
!  Get the job type
!  ----------------
   jobType = getTaskType(meshName, no_of_solutions, solutionNames, solutionTypes, fixedOrder, Nout, basis)
!
!  Construct Spectral basis
!  ------------------------
   call ConstructSpectralBasis
!
!  Perform the conversion to tecplot
!  ---------------------------------
   select case (jobType)
   case (MESH_2_PLT)
      call Mesh2Plt(meshName)

   case (SOLUTION_2_PLT)
      do iSol = 1, no_of_solutions

         select case (solutionTypes(iSol))
         case ( SOLUTION_FILE )
            call Solution2Plt(meshName, solutionNames(iSol), fixedOrder, basis, Nout)        

         case ( SOLUTION_AND_GRADIENTS_FILE )
            call Solution2Plt(meshName, solutionNames(iSol), fixedOrder, basis, Nout)        

         case ( STATS_FILE )
            call Stats2Plt(meshName, solutionNames(iSol), fixedOrder, basis, Nout)

         end select
      end do

   case (UNKNOWN_JOB)
      call exit(UNKNOWN_JOB)

   end select

end program horses2plt
