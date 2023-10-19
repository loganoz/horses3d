!
!/////////////////////////////////////////////////////////////////////////////////////
!
!     This subroutines handle the creation and reading of HORSES unformatted files.
!  There are four types of files:
!        * Mesh files: stores coordinates of interpolation points (for the initial mesh)
!        * Solution files: stores the conservative variables Q for all elements
!        * Solution and gradient files: stores the conservative variables and their gradients.
!        * Statistics files: stores the statistics values.
!        * Zone solution files: stores the conservative variables Q and its temporal derivate for the faces of a zone
!
!  These files always start as follows:
!
!        1- character variable (128) containing the name of the file
!        2- integer variable containing the file type
!        3- integer containing the node type (GAUSS/GAUSS-LOBATTO)
!        4- integer variable containing the number of elements
!        5- the iteration
!        6- the time
!        7- reference quantities (6 real numbers for gamma, R, V, Rho, T, Mach)
!        8- a integer termination indicator (BEGINNING_DATA)
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
#ifdef _HAS_MPI_
   use mpi
#endif
   
   private
   public      :: MESH_FILE, SOLUTION_FILE, SOLUTION_AND_GRADIENTS_FILE, STATS_FILE, ZONE_MESH_FILE
   public      :: ZONE_SOLUTION_FILE, ZONE_SOLUTION_AND_DOT_FILE
   public      :: SOLUTION_AND_SENSOR_FILE, SOLUTION_AND_GRADIENTS_AND_SENSOR_FILE
   public      :: BEGINNING_DATA
   public      :: SOLFILE_STR_LEN, POS_INIT_DATA
   public      :: NO_OF_SAVED_REFS, GAMMA_REF, RGAS_REF, V_REF, RHO_REF, T_REF, MACH_REF, RE_REF

   public      :: CreateNewSolutionFile, writeArray, SealSolutionFile, getSolutionFileType
   public      :: getSolutionFileNoOfElements, getSolutionFileName
   public      :: getSolutionFileDataInitFlag
   public      :: getSolutionFileArrayDimensions, getSolutionFileReferenceValues
   public      :: getSolutionFileNodeType, getSolutionFileTimeAndIteration
   public      :: putSolutionFileInReadDataMode, putSolutionFileInWriteDataMode
   
   public      :: POS_FILETYPE
!
!  Possible solution file types
!  ----------------------------
   integer, parameter      :: MESH_FILE                              = 1
   integer, parameter      :: SOLUTION_FILE                          = 2
   integer, parameter      :: SOLUTION_AND_GRADIENTS_FILE            = 3
   integer, parameter      :: STATS_FILE                             = 4
   integer, parameter      :: ZONE_MESH_FILE                         = 5
   integer, parameter      :: ZONE_SOLUTION_FILE                     = 6
   integer, parameter      :: ZONE_SOLUTION_AND_DOT_FILE             = 7
   integer, parameter      :: SOLUTION_AND_SENSOR_FILE               = 8
   integer, parameter      :: SOLUTION_AND_GRADIENTS_AND_SENSOR_FILE = 9

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
!
!  ------------
!  File offsets
!  ------------
!
   integer, parameter   :: POS_FILENAME     = 1
   integer, parameter   :: POS_FILETYPE     = POS_FILENAME + SOLFILE_STR_LEN * SIZEOF_CHAR
   integer, parameter   :: POS_NODETYPE     = POS_FILETYPE + SIZEOF_INT
   integer, parameter   :: POS_NOOFELEMENTS = POS_NODETYPE + SIZEOF_INT
   integer, parameter   :: POS_ITER         = POS_NOOFELEMENTS + SIZEOF_INT
   integer, parameter   :: POS_TIME         = POS_ITER + SIZEOF_INT
   integer, parameter   :: POS_REFS         = POS_TIME + SIZEOF_RP
   integer, parameter   :: POS_TERMINATOR   = POS_REFS + NO_OF_SAVED_REFS * SIZEOF_RP
   integer, parameter   :: POS_INIT_DATA    = POS_TERMINATOR + SIZEOF_INT

   interface writeArray
      module procedure  :: write0DArray, write1DArray, write2DArray, write3DArray
      module procedure  :: write4DArray, write5DArray
   end interface writeArray
    
   contains
      subroutine CreateNewSolutionFile(name, type_, nodes, no_of_elements, iter, time, refs)
!
!        **********************************************************************
!           This function creates a new solution file and writes its header
!        **********************************************************************
!
         use MPI_Process_Info
         implicit none
         character(len=*), intent(in)              :: name
         integer,          intent(in)              :: type_
         integer,          intent(in)              :: nodes
         integer,          intent(in)              :: no_of_elements
         integer,          intent(in)              :: iter
         real(kind=RP),    intent(in)              :: time
         real(kind=RP),    intent(in)              :: refs(NO_OF_SAVED_REFS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                          :: fID, ierr
         character(len=SOLFILE_STR_LEN)   :: auxstr
!
!        Only root creates the file
!        --------------------------
         if ( MPI_Process % isRoot ) then
!
!           Check consistency
!           -----------------         
            select case (type_)
            case(MESH_FILE)
            case(SOLUTION_FILE)
            case(SOLUTION_AND_GRADIENTS_FILE)
            case(STATS_FILE)
            case(ZONE_MESH_FILE)
            case(ZONE_SOLUTION_FILE)
            case(ZONE_SOLUTION_AND_DOT_FILE)
            case(SOLUTION_AND_SENSOR_FILE)
            case(SOLUTION_AND_GRADIENTS_AND_SENSOR_FILE)
            case default
               print*, "Incorrect solution file type", type_
               errorMessage(STD_OUT)
               error stop
            end select
!   
!           Open the file
!           -------------
            open(newunit=fID, file=trim(name), action="write", &
                 status="replace", form="unformatted", access = "stream") 
!   
!           Write the title
!           ---------------
            auxstr = name
            write(fID, POS=POS_FILENAME) trim(auxstr)
!   
!           Write the type
!           --------------
            write(fID, POS=POS_FILETYPE) type_         
!   
!           Write the nodes type
!           --------------------
            write(fID, POS=POS_NODETYPE) nodes
!   
!           Write the number of elements
!           ----------------------------
            write(fID, POS=POS_NOOFELEMENTS) no_of_elements
!   
!           Save the current time and iteration
!           -----------------------------------
            write(fID, POS=POS_ITER) iter
            write(fID, POS=POS_TIME) time
!   
!           Save refernence and gas data
!           ----------------------------
            write(fID, POS=POS_REFS) refs
!   
!           Introduce the terminator indicator
!           ----------------------------------                    ! TODO: Write which physics are being used
            write(fID, POS=POS_TERMINATOR) BEGINNING_DATA
!   
!           Close the file
!           --------------
            close(fid)

         end if
!
!        Introduce a barrier
!        -------------------
#ifdef _HAS_MPI_
         call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

      end subroutine CreateNewSolutionFile

      subroutine SealSolutionFile(name)
!
!        **********************************************
!        This function appends a terminator to the file
!        to mark when the data has finished
!        **********************************************
!
         use MPI_Process_Info
         implicit none
         character(len=*), intent(in)              :: name
         integer                                   :: fid, ierr
         integer(kind=AddrInt)                     :: pos
!
!        Add a barrier to make sure that all processes have written their data
!        ---------------------------------------------------------------------
#ifdef _HAS_MPI_
         call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

         if ( MPI_Process % isRoot ) then
            open(newunit=fID, file=trim(name), action="write", status="old", &
                 form="unformatted", access = "stream") 
            inquire(unit=fid, size=pos) 
            write(fID, pos=pos+1_AddrInt) END_OF_FILE
            close(fID)
         end if

      end subroutine SealSolutionFile

      character(len=SOLFILE_STR_LEN) function getSolutionFileName(fileName)
         implicit none
         character(len=*), intent(in)  :: fileName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: fid
!
!        Open file
!        ---------
         open(newunit=fid, file=trim(fileName), status="old", action="read", form="unformatted", access="stream")
!
!        Get the file type
!        -----------------
         read(fid, pos=POS_FILENAME) getSolutionFileName
!
!        Close file
!        ----------
         close(fid)

      end function getSolutionFileName
   
      integer function getSolutionFileType(fileName)
         implicit none
         character(len=*), intent(in)  :: fileName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: fid
         character(len=128)   :: cosa
         integer              :: cosa2
!
!        Open file
!        ---------
         open(newunit=fid, file=trim(fileName), status="old", action="read", form="unformatted", access="stream")
!
!        Get the file type
!        -----------------
         read(fid, pos=POS_FILETYPE) getSolutionFileType
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
         integer  :: fid
!
!        Open file
!        ---------
         open(newunit=fid, file=trim(fileName), status="old", action="read", form="unformatted", access="stream")
!
!        Get the number of elements
!        --------------------------
         read(fid, pos=POS_NOOFELEMENTS) getSolutionFileNoOfElements
!
!        Close file
!        ----------
         close(fid)

      end function getSolutionFileNoOfElements

      integer function getSolutionFileNodeType(fileName)
         implicit none
         character(len=*), intent(in)  :: fileName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: fid
!
!        Open file
!        ---------
         open(newunit=fid, file=trim(fileName), status="old", action="read", form="unformatted", access="stream")
!
!        Get the node type
!        -----------------
         read(fid, pos=POS_NODETYPE) getSolutionFileNodeType
!
!        Close file
!        ----------
         close(fid)

      end function getSolutionFileNodeType

      subroutine getSolutionFileTimeAndIteration(fileName, iter, time)
         implicit none
         character(len=*), intent(in)  :: fileName
         integer,          intent(out) :: iter
         real(kind=RP),    intent(out) :: time
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: fid
!
!        Open file
!        ---------
         open(newunit=fid, file=trim(fileName), status="old", action="read", form="unformatted", access="stream")
!  
!        Get iteration
!        -------------
         read(fid, pos=POS_ITER) iter
!
!        Get time
!        --------
         read(fid, pos=POS_TIME) time
!
!        Close file
!        ----------
         close(fid)

      end subroutine getSolutionFileTimeAndIteration

      integer function getSolutionFileDataInitFlag(fileName)
         implicit none
         character(len=*), intent(in)  :: fileName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: fid
!
!        Open file
!        ---------
         open(newunit=fid, file=trim(fileName), status="old", action="read", form="unformatted", access="stream")
!
!        Get terminator
!        --------------
         read(fid, pos=POS_TERMINATOR) getSolutionFileDataInitFlag
!
!        Close file
!        ----------
         close(fid)

      end function getSolutionFileDataInitFlag

      function getSolutionFileReferenceValues(fileName)
         implicit none
         character(len=*), intent(in)               :: fileName
         real(kind=RP), dimension(NO_OF_SAVED_REFS) :: getSolutionFileReferenceValues
!
!        ---------------
!        Local variables
!        ---------------
!
         integer :: fid
!
!        Open file
!        ---------
         open(newunit=fid, file=trim(fileName), status="old", action="read", form="unformatted", access="stream")
!
!        Get dimensionless and reference values
!        --------------------------------------
         read(fid, pos=POS_REFS) getSolutionFileReferenceValues
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
         integer :: fid, flag
!
!        Open file
!        ---------
         open(newunit=fid, file=trim(fileName), status="old", action="read", &
                           form="unformatted" , access="stream")

!
!        Move to the INIT_DATA position
!        ------------------------------
         read(fid, pos=POS_TERMINATOR) flag
      
         if ( flag .ne. BEGINNING_DATA ) then
            print*, "Wrong beginning data specifier"
            errorMessage(STD_OUT)
            error stop
         end if
!
!        Return the file ID
!        ------------------
         putSolutionFileInReadDataMode = fid

      end function putSolutionFileInReadDataMode

      integer function putSolutionFileInWriteDataMode(fileName)
         implicit none
         character(len=*), intent(in)  :: fileName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer :: fid, flag
!
!        Open file
!        ---------
         open(newunit=fid, file=trim(fileName), status="old", action="readwrite", &
                           form="unformatted" , access="stream")

!
!        Move to the INIT_DATA position
!        ------------------------------
         read(fid, pos=POS_TERMINATOR) flag

         if ( flag .ne. BEGINNING_DATA ) then
            print*, "Wrong beginning data specifier"
            errorMessage(STD_OUT)
            error stop
         end if
!
!        Return the file ID
!        ------------------
         putSolutionFileInWriteDataMode = fid

      end function putSolutionFileInWriteDataMode
!
!/////////////////////////////////////////////////////////////////////
!
!     Write array procedures
!     ----------------------
!
!/////////////////////////////////////////////////////////////////////
!
      subroutine Write0DArray(fID, array, position)
         implicit none
         integer, intent(in)                         :: fID
         real(kind=RP), intent(in)                   :: array
         integer(kind=AddrInt), intent(in), optional :: position

         if ( present(position) ) then
            write(fID, pos = position) 0
            write(fID) 1
            write(fID) array

         else
            write(fID) 0
            write(fID) 1
            write(fID) array

         end if

      end subroutine Write0DArray

      subroutine Write1DArray(fID,array, position)
         implicit none
         integer, intent(in)                         :: fID
         real(kind=RP), intent(in)                   :: array(:)
         integer(kind=AddrInt), intent(in), optional :: position

         if ( present(position) ) then
            write(fID) 1
            write(fID) size(array,1)
            write(fID) array
   
         else
            write(fID) 1
            write(fID) size(array,1)
            write(fID) array
         end if

      end subroutine Write1DArray

      subroutine Write2DArray(fID,array, position)
         implicit none
         integer, intent(in)                         :: fID
         real(kind=RP), intent(in)                   :: array(:,:)
         integer(kind=AddrInt), intent(in), optional :: position

         if ( present(position) ) then
            write(fID, pos = position) 2
            write(fID) size(array,1), size(array,2)
            write(fID) array
   
         else
            write(fID) 2
            write(fID) size(array,1), size(array,2)
            write(fID) array
         end if

      end subroutine Write2DArray

      subroutine Write3DArray(fID,array, position)
         implicit none
         integer, intent(in)                         :: fID
         real(kind=RP), intent(in)                   :: array(:,:,:)
         integer(kind=AddrInt), intent(in), optional :: position

         if ( present(position) ) then
            write(fID, pos = position) 3
            write(fID) size(array,1), size(array,2), size(array,3)
            write(fID) array
   
         else
            write(fID) 3
            write(fID) size(array,1), size(array,2), size(array,3)
            write(fID) array
         end if

      end subroutine Write3DArray

      subroutine Write4DArray(fID,array, position)
         implicit none
         integer, intent(in)                         :: fID
         real(kind=RP), intent(in)                   :: array(:,:,:,:)
         integer(kind=AddrInt), intent(in), optional :: position

         if ( present(position) ) then
            write(fID, pos = position) 4
            write(fID) size(array,1), size(array,2), size(array,3), size(array,4)
            write(fID) array
   
         else
            write(fID) 4
            write(fID) size(array,1), size(array,2), size(array,3), size(array,4)
            write(fID) array
         end if

      end subroutine Write4DArray

      subroutine Write5DArray(fID,array, position)
         implicit none
         integer, intent(in)                         :: fID
         real(kind=RP), intent(in)                   :: array(:,:,:,:,:)
         integer(kind=AddrInt), intent(in), optional :: position

         if ( present(position) ) then
            write(fID, pos = position) 5
            write(fID) size(array,1), size(array,2), size(array,3), size(array,4), size(array,5)
            write(fID) array
   
         else
            write(fID) 5
            write(fID) size(array,1), size(array,2), size(array,3), size(array,4), size(array,5)
            write(fID) array
         end if

      end subroutine Write5DArray
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     Read array procedures
!     ---------------------
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine getSolutionFileArrayDimensions(fid,N,pos)
         implicit none
         integer, intent(in)           :: fid
         integer, intent(out)          :: N(:)
         integer, intent(in), optional :: pos
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: arrayDimension

         if ( present(pos) ) then
            read(fid, pos=pos) arrayDimension
   
         else
            read(fid) arrayDimension

         end if

         if ( size(N) .ne. arrayDimension) then
            print*, "Array found in file dimensions does not match that of the introduced variable. File dimensions: ", &
                    arrayDimension, ", Variable: ", size(N)
            errorMessage(STD_OUT)
            error stop
         end if

         read(fid) N

      end subroutine getSolutionFileArrayDimensions

end module SolutionFile