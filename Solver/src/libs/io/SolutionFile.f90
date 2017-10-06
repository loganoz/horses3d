!
!/////////////////////////////////////////////////////////////////////////////////////
!
!     This subroutines handle the creation and reading of HORSES unformatted files.
!  There are four types of files:
!        * Mesh files: stores coordinates of interpolation points (for the initial mesh)
!        * Solution files: stores the conservative variables Q for all elements
!        * Solution and gradient files: stores the conservative variables and their gradients.
!        * Statistics files: stores the statistics values.
!
!  These files always start as follows:
!
!        1- character variable (128) containing the name of the file
!        2- integer variable containing the file type
!        3- integer variable containing the number of elements
!        4- (solution files only) the iteration
!        5- (solution files only) the time
!        6- a integer termination indicator (BEGINNING_DATA)
!
!  Next, the data is stored. This is always stored in the way:
!
!     ** for each element:
!        1- Integer containing the number of dimensions
!        2- Integers containing array dimensions: (e.g.) Nx,Ny,Nz,NCONS
!        3- Real array containing the data: (e.g.) Q
!
!  Last, the END_OF_FILE indicator (integer) is introduced (just to ensure)
!
!/////////////////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
module SolutionFiles
   use SMConstants
   
   private
   public      :: MESH_FILE, SOLUTION_FILE, SOLUTION_AND_GRADIENTS_FILE, STATS_FILE
   public      :: SOLFILE_STR_LEN

   public      :: CreateNewSolutionFile, writeArray
!
!  Possible solution file types
!  ----------------------------
   integer, parameter      :: MESH_FILE = 1
   integer, parameter      :: SOLUTION_FILE = 2
   integer, parameter      :: SOLUTION_AND_GRADIENTS_FILE = 3
   integer, parameter      :: STATS_FILE = 4  

   integer, parameter      :: SOLFILE_STR_LEN = 128
   integer, parameter      :: END_OF_FILE    = 99
   integer, parameter      :: BEGINNING_DATA = 88

   interface writeArray
      module procedure  :: write0DArray, write1DArray, write2DArray, write3DArray
      module procedure  :: write4DArray, write5DArray
   end interface writeArray
    
   contains
      integer function CreateNewSolutionFile(name, type_, no_of_elements, iter, time)
         implicit none
         character(len=*), intent(in)              :: name
         integer,          intent(in)              :: type_
         integer,          intent(in)              :: no_of_elements
         integer,          intent(in), optional    :: iter
         real(kind=RP),    intent(in), optional    :: time
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                          :: fID
         character(len=SOLFILE_STR_LEN)   :: auxstr
!
!        Check consistency
!        -----------------         
         select case (type_)
         case(MESH_FILE)
         case(SOLUTION_FILE)
         case(SOLUTION_AND_GRADIENTS_FILE)
         case(STATS_FILE)
         case default
            print*, "Incorrect solution file type"
            errorMessage(STD_OUT)
            stop
         end select
!
!        Open the file
!        -------------
         open(newunit=fID, file=trim(name), action="write", status="unknown", form="unformatted") 
!
!        Write the title
!        ---------------
         auxstr = name
         write(fID) auxstr
!
!        Write the type
!        --------------
         write(fID) type_         
!
!        Write the number of elements
!        ----------------------------
         write(fID) no_of_elements
!
!        If the file is solution file, save the current time and iteration
!        -----------------------------------------------------------------
         if ( type_ .ne. MESH_FILE ) then
            if ( (.not. present(time)) .or. (.not. present(iter)) ) then
               print*, "Missing time and/or iteration values"
               errorMessage(STD_OUT)
               stop
            end if
            write(fID) iter
            write(fID) time
         end if
!
!        Introduce the terminator indicator
!        ----------------------------------
         write(fID) BEGINNING_DATA
!
!        Return the file ID
!        ------------------         
         CreateNewSolutionFile = fID

      end function CreateNewSolutionFile

      subroutine CloseSolutionFile(fID)
         implicit none
         integer, intent(in)     :: fID
   
         write(fID) END_OF_FILE
         close(fID)

      end subroutine CloseSolutionFile

      subroutine Write0DArray(fID, array)
         implicit none
         integer, intent(in)           :: fID
         real(kind=RP), intent(in)     :: array

         write(fID) 0
         write(fID) 1
         write(fID) array

      end subroutine Write0DArray

      subroutine Write1DArray(fID,array)
         implicit none
         integer, intent(in)           :: fID
         real(kind=RP), intent(in)     :: array(:)

         write(fID) 1
         write(fID) size(array,1)
         write(fID) array

      end subroutine Write1DArray

      subroutine Write2DArray(fID,array)
         implicit none
         integer, intent(in)           :: fID
         real(kind=RP), intent(in)     :: array(:,:)

         write(fID) 2
         write(fID) size(array,1), size(array,2)
         write(fID) array

      end subroutine Write2DArray

      subroutine Write3DArray(fID,array)
         implicit none
         integer, intent(in)           :: fID
         real(kind=RP), intent(in)     :: array(:,:,:)

         write(fID) 3
         write(fID) size(array,1), size(array,2), size(array,3)
         write(fID) array

      end subroutine Write3DArray

      subroutine Write4DArray(fID,array)
         implicit none
         integer, intent(in)           :: fID
         real(kind=RP), intent(in)     :: array(:,:,:,:)

         write(fID) 4
         write(fID) size(array,1), size(array,2), size(array,3), size(array,4)
         write(fID) array

      end subroutine Write4DArray

      subroutine Write5DArray(fID,array)
         implicit none
         integer, intent(in)           :: fID
         real(kind=RP), intent(in)     :: array(:,:,:,:,:)

         write(fID) 5
         write(fID) size(array,1), size(array,2), size(array,3), size(array,4), size(array,5)
         write(fID) array

      end subroutine Write5DArray



end module SolutionFiles
