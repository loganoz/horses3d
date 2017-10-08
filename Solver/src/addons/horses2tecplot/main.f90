program horses2plt
   use SMConstants
   use getTask
   use Mesh2PltModule
   implicit none
   integer                    :: jobType
   character(len=LINE_LENGTH) :: meshName
   character(len=LINE_LENGTH) :: solutionName
!
!  Get the job type
!  ----------------
   jobType = getTaskType(meshName,solutionName)
!
!  Perform the conversion to tecplot
!  ---------------------------------
   select case (jobType)
   case (MESH_2_PLT)
      call Mesh2Plt(meshName)

   case (UNKNOWN_JOB)
      call exit(UNKNOWN_JOB)

   end select










end program horses2plt
