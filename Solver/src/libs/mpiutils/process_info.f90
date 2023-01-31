module MPI_Process_Info
   use SMConstants
#ifdef _HAS_MPI_
   use mpi
#endif

   private
   public   MPI_Process, DEFAULT_TAG

   integer, parameter   :: DEFAULT_TAG = 99

   type MPI_Process_t
      logical     :: doMPIRootAction = .false.
      logical     :: doMPIAction     = .false.
      logical     :: isRoot          = .false.
      integer     :: nProcs          = 1
      integer     :: rank            = 0
      contains
         procedure   :: Init => MPI_Process_Init
         procedure   :: Close => MPI_Process_Close
   end type MPI_Process_t

   type(MPI_Process_t)  :: MPI_Process

   contains
      subroutine MPI_Process_Init(self)
         implicit none
         class(MPI_Process_t),   intent(out)  :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: ierr

#ifdef _HAS_MPI_   
         call mpi_init(ierr)
         call mpi_comm_rank(MPI_COMM_WORLD, self % rank, ierr)
         call mpi_comm_size(MPI_COMM_WORLD, self % nProcs, ierr)
#else
         self % rank = 0
         self % nProcs = 1
#endif
!
!        ----------------------------------
!        The rank 0 will be defined as root
!        ----------------------------------
!
         if ( self % rank .eq. 0 ) self % isRoot = .true.
!
!        ------------------------------------------------------------
!        doMPIAction is just true if the code is running with MPI
!        ------------------------------------------------------------
!
         if ( self % nProcs .ne. 1 ) self % doMPIAction = .true.
!
!        ------------------------------------------------------------
!        The doMPIRootAction is a variable to control root operations
!        that are NOT required without MPI (or with 1 process)
!        ------------------------------------------------------------
!
         if ( (self % nProcs .ne. 1) .and. (self % isRoot) ) then
            self % doMPIRootAction = .true.
         end if

      end subroutine MPI_Process_Init

      subroutine MPI_Process_Close(self)
         implicit none
         class(MPI_Process_t),   intent(in)  :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: ierr

#ifdef _HAS_MPI_
         call mpi_finalize(ierr)
#endif 
      end subroutine MPI_Process_Close 

end module MPI_Process_Info