!
!/////////////////////////////////////////////////////////////////////
!
!     This module explores the command line flags introduced
!  by the user. Software usage is:
!
!     horses2tecplot *.hmesh *.hsol --interpolate=Npoints
!
!  where the flag --interpolate is optional. If it is not present,
!  the solution will be exported using the original Gauss points.
!
!/////////////////////////////////////////////////////////////////////
!
module getTask
   use SMConstants
   implicit none
   
   private
   public   MESH_2_PLT, SOLUTION_2_PLT, UNKNOWN_JOB

   public   getTaskType

   integer, parameter   :: UNKNOWN_JOB = 0
   integer, parameter   :: MESH_2_PLT = 1
   integer, parameter   :: SOLUTION_2_PLT = 2

   contains

      integer function getTaskType(meshName, no_of_solutions, solutionNames, solutionTypes, performInterpolation, Npoints)
         use SolutionFile
         implicit none
         character(len=*),                        intent(out) :: meshName
         integer,                                 intent(out) :: no_of_solutions
         character(len=LINE_LENGTH), allocatable, intent(out) :: solutionNames(:)
         integer, allocatable,                    intent(out) :: solutionTypes(:)
         logical,                                 intent(out) :: performInterpolation
         integer,                                 intent(out) :: Npoints
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: no_of_arguments
         integer     :: i, sol
         logical     :: meshFilePresent
         character(len=LINE_LENGTH) :: auxiliarName
!
!        Get number of command arguments         
!        -------------------------------
         no_of_arguments = command_argument_count()
!
!        Exit if no input arguments are specified
!        ----------------------------------------
         if ( no_of_arguments .eq. 0 ) then
            write(STD_OUT,'(A)') "No mesh file and/or solution file specified"
            getTaskType = UNKNOWN_JOB
            return
         end if
!
!        Check if the mesh file is present
!        ---------------------------------
         meshFilePresent = .false.
         do i = 1, no_of_arguments
            call get_command_argument(i, meshName)
            if ( getSolutionFileType(meshName) .eq. MESH_FILE ) then
               meshFilePresent = .true.
               exit
            end if 
         end do
!
!        Exit if the mesh file is not present
!        ------------------------------------
         if ( .not. meshFilePresent ) then
            write(STD_OUT,'(A)') "Mesh file was not specified"
            getTaskType = UNKNOWN_JOB
            return
         end if
!
!        Check if the solution file is present
!        -------------------------------------
         no_of_solutions = 0
         do i = 1, no_of_arguments
            call get_command_argument(i, auxiliarName)
            if ( getSolutionFileType(auxiliarName) .ne. MESH_FILE ) then
               no_of_solutions = no_of_solutions + 1 
            end if
         end do

         if ( no_of_solutions .ne. 0 ) then
            allocate( solutionNames(no_of_solutions) )
            allocate( solutionTypes(no_of_solutions) )
            
            sol = 0
            do i = 1, no_of_arguments
               call get_command_argument(i, auxiliarName)
               if ( getSolutionFileType(auxiliarName) .ne. MESH_FILE ) then   
                  sol = sol + 1 
                  solutionNames(sol) = trim(auxiliarName)
                  solutionTypes(sol) = getSolutionFileType(auxiliarName)
               end if
            end do

         end if
!
!        Select the job type
!        -------------------             
         if ( no_of_solutions .ne. 0 ) then
            getTaskType = SOLUTION_2_PLT
   
         else
            getTaskType = MESH_2_PLT

         end if

         performInterpolation = .false.
         Npoints = 0
         
      end function getTaskType
end module getTask
