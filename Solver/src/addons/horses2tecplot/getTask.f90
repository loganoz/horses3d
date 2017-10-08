module getTask
   use SMConstants
   implicit none
   
   private
   public   MESH_2_PLT, SOLUTION_2_PLT, SOLUTION_AND_GRADIENTS_2_PLT
   public   STATS_2_PLT, UNKNOWN_JOB

   public   getTaskType

   integer, parameter   :: UNKNOWN_JOB = 0
   integer, parameter   :: MESH_2_PLT = 1
   integer, parameter   :: SOLUTION_2_PLT = 2
   integer, parameter   :: SOLUTION_AND_GRADIENTS_2_PLT = 3
   integer, parameter   :: STATS_2_PLT = 4

   contains

      integer function getTaskType(meshName, solutionName)
         use SolutionFile
         implicit none
         character(len=*), intent(out)    :: meshName
         character(len=*), intent(out)    :: solutionName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: no_of_arguments

         no_of_arguments = command_argument_count()

         select case ( no_of_arguments )
         case (0)
            write(STD_OUT,'(A)') "No mesh file and solution file selected" 
            getTaskType = UNKNOWN_JOB

         case (1)
            call get_command_argument(1, meshName) 
            
            if ( getSolutionFileType(meshName) .eq. MESH_FILE ) then
               getTaskType = MESH_2_PLT
            else
               write(STD_OUT,'(A,A,A)') 'File "',trim(meshName),'" is not a mesh file'
               getTaskType = UNKNOWN_JOB
            end if
         end select
         
      end function getTaskType

end module getTask
