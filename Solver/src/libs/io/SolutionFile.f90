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
module SolutionFile
   use SMConstants
   
   private
   public      :: MESH_FILE, SOLUTION_FILE, SOLUTION_AND_GRADIENTS_FILE, STATS_FILE
   public      :: BEGINNING_DATA
   public      :: SOLFILE_STR_LEN
   public      :: NO_OF_SAVED_REFS, GAMMA_REF, RGAS_REF, V_REF, RHO_REF, T_REF, MACH_REF

   public      :: CreateNewSolutionFile, writeArray, CloseSolutionFile, getSolutionFileType
   public      :: putSolutionFileInReadDataMode, getSolutionFileNoOfElements
   public      :: getSolutionFileArrayDimensions, getSolutionFileReferenceValues
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

   integer, parameter      :: NO_OF_SAVED_REFS = 6
   integer, parameter      :: GAMMA_REF = 1
   integer, parameter      :: RGAS_REF  = 2
   integer, parameter      :: V_REF     = 3
   integer, parameter      :: RHO_REF   = 4
   integer, parameter      :: T_REF     = 5
   integer, parameter      :: MACH_REF  = 6

   interface writeArray
      module procedure  :: write0DArray, write1DArray, write2DArray, write3DArray
      module procedure  :: write4DArray, write5DArray
   end interface writeArray
    
   contains
      integer function CreateNewSolutionFile(name, type_, no_of_elements, iter, time, refs)
         implicit none
         character(len=*), intent(in)              :: name
         integer,          intent(in)              :: type_
         integer,          intent(in)              :: no_of_elements
         integer,          intent(in)              :: iter
         real(kind=RP),    intent(in)              :: time
         real(kind=RP),    intent(in)              :: refs(NO_OF_SAVED_REFS)
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
!        Save the current time and iteration
!        -----------------------------------
         write(fID) iter
         write(fID) time
!
!        Save refernence and gas data
!        ----------------------------
         write(fID) refs
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
   
      integer function getSolutionFileType(fileName)
         implicit none
         character(len=*), intent(in)  :: fileName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: fid
         character(len=SOLFILE_STR_LEN)   :: fileNameInFile
!
!        Open file
!        ---------
         open(newunit=fid, file=trim(fileName), status="old", action="read", form="unformatted")
!
!        Get the file name
!        -----------------
         read(fid) fileNameInFile
!
!        Get the file type
!        -----------------
         read(fid) getSolutionFileType
!
!        Close file
!        ----------
         close(fid)

      end function getSolutionFileType

      integer function getSolutionFileNoOfElements(fileName)
         implicit none
         character(len=*), intent(in)  :: fileName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: fid, no_of_elements, initial_iter, fileType, flag
         real(kind=RP)  :: initial_time
         character(len=SOLFILE_STR_LEN)   :: fileNameInFile
!
!        Open file
!        ---------
         open(newunit=fid, file=trim(fileName), status="old", action="read", form="unformatted")
!
!        Get the file name
!        -----------------
         read(fid) fileNameInFile
!
!        Get the file type
!        -----------------
         read(fid) fileType
!
!        Get the number of elements
!        --------------------------
         read(fid) getSolutionFileNoOfElements
!
!        Close file
!        ----------
         close(fid)

      end function getSolutionFileNoOfElements

      function getSolutionFileReferenceValues(fileName)
         implicit none
         character(len=*), intent(in)               :: fileName
         real(kind=RP), dimension(NO_OF_SAVED_REFS) :: getSolutionFileReferenceValues
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: fid, no_of_elements, initial_iter, fileType, flag
         real(kind=RP)  :: initial_time
         character(len=SOLFILE_STR_LEN)   :: fileNameInFile
!
!        Open file
!        ---------
         open(newunit=fid, file=trim(fileName), status="old", action="read", form="unformatted")
!
!        Get the file name
!        -----------------
         read(fid) fileNameInFile
!
!        Get the file type
!        -----------------
         read(fid) fileType
!
!        Get the number of elements
!        --------------------------
         read(fid) 
!
!        Skip the time and iteration
!        ---------------------------
         read(fid)
         read(fid)
!
!        Get dimensionless and reference values
!        --------------------------------------
         read(fid) getSolutionFileReferenceValues
!
!        Close file
!        ----------
         close(fid)

      end function getSolutionFileReferenceValues

      integer function putSolutionFileInReadDataMode(fileName)
         implicit none
         character(len=*), intent(in)  :: fileName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: fid, no_of_elements, initial_iter, fileType, flag
         real(kind=RP)  :: initial_time
         character(len=SOLFILE_STR_LEN)   :: fileNameInFile
!
!        Open file
!        ---------
         open(newunit=fid, file=trim(fileName), status="old", action="read", form="unformatted")
!
!        Get the file name
!        -----------------
         read(fid) fileNameInFile
!
!        Get the file type
!        -----------------
         read(fid) fileType
!
!        Get the number of elements
!        --------------------------
         read(fid) no_of_elements
!
!        Get the initial time and iteration
!        ----------------------------------
         read(fid) initial_iter
         read(fid) initial_time        
!
!        Get the reference values
!        ------------------------
         read(fid)
!
!        Get the beginning data mark
!        ---------------------------
         read(fid) flag

         if ( flag .ne. BEGINNING_DATA ) then
            write(STD_OUT,'(A,A,A)') 'Error reading file "',trim(fileName),'". The beginning data mark was not found.'
            close(fid)
            putSolutionFileInReadDataMode = 0
         else
            putSolutionFileInReadDataMode = fid
         end if

      end function putSolutionFileInReadDataMode
!
!/////////////////////////////////////////////////////////////////////
!
!     Write array procedures
!     ----------------------
!
!/////////////////////////////////////////////////////////////////////
!
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
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     Read array procedures
!     ---------------------
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine getSolutionFileArrayDimensions(fid,N)
         implicit none
         integer, intent(in)     :: fid
         integer, intent(out)    :: N(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: arrayDimension

         read(fid) arrayDimension

         if ( size(N) .ne. arrayDimension) then
            print*, "Array found in file dimensions does not match that of the introduced variable."
            errorMessage(STD_OUT)
            stop
         end if

         read(fid) N

      end subroutine getSolutionFileArrayDimensions

end module SolutionFile
