#include "Includes.h"
module MPI_Face_Class
   use SMConstants
   use MPI_Process_Info
#ifdef _HAS_MPI_ 
   use mpi 
#endif 
   implicit none

   private
   public  MPI_Face_t, mpi_faces, MPI_Faces_Constructed

   public  ConstructMPIFaces, DestructMPIFaces
   public  ConstructMPIFacesStorage

   type MPI_Face_t
      integer                    :: no_of_faces
      integer, allocatable       :: faceIDs(:)
      integer, allocatable       :: elementSide(:)
      integer                    :: Qrecv_req
      integer                    :: gradQrecv_req
      integer                    :: sizeQ
      integer                    :: sizeU_xyz
      real(kind=RP), allocatable :: Qsend(:)
      real(kind=RP), allocatable :: U_xyzsend(:)
      real(kind=RP), allocatable :: Qrecv(:)
      real(kind=RP), allocatable :: U_xyzrecv(:)
      contains
         procedure   :: Construct        => MPI_Face_Construct
         procedure   :: Destruct         => MPI_Face_Destruct
         procedure   :: SendQ            => MPI_Face_SendQ
         procedure   :: RecvQ            => MPI_Face_RecvQ
         procedure   :: SendU_xyz        => MPI_Face_SendU_xyz
         procedure   :: RecvU_xyz        => MPI_Face_RecvU_xyz
         procedure   :: WaitForSolution  => MPI_Face_WaitForSolution
         procedure   :: WaitForGradients => MPI_Face_WaitForGradients
   end type MPI_Face_t

   logical                          :: MPI_Faces_Constructed = .false.
   type(MPI_Face_t), allocatable    :: mpi_faces(:)

   interface MPI_Face_t
      module procedure Construct_MPI_Face
   end interface MPI_Face_t

   contains

      subroutine ConstructMPIFaces()
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: domain

         if ( MPI_Process % doMPIAction ) then
            allocate(mpi_faces(MPI_Process % nProcs))

            do domain = 1, MPI_Process % nProcs
               mpi_faces(domain) = MPI_Face_t()
            end do
         end if

         MPI_Faces_Constructed = .true.

      end subroutine ConstructMPIFaces

      subroutine ConstructMPIFacesStorage(NCONS, NGRAD, NDOFS)
!
!        ***************************************************
!           Allocates buffers to send and receive.
!           This subroutine is called from:
!               mesh % SetElementConnectivitiesAndLinkFaces
!        ***************************************************
!
         implicit none
         integer, intent(in)     :: NCONS, NGRAD
         integer, intent(in)     :: NDOFS(MPI_Process % nProcs)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: domain

         do domain = 1, MPI_Process % nProcs
            mpi_faces(domain) % sizeQ     = NCONS * NDOFS(domain)
            mpi_faces(domain) % sizeU_xyz = 3 * NGRAD * NDOFS(domain)

            if ( NDOFS(domain) .gt. 0 ) then
               allocate( mpi_faces(domain) % Qsend(NCONS * NDOFS(domain)) )
               allocate( mpi_faces(domain) % U_xyzsend(3 * NGRAD * NDOFS(domain)) )
               allocate( mpi_faces(domain) % Qrecv(NCONS * NDOFS(domain)) )
               allocate( mpi_faces(domain) % U_xyzrecv(3 * NGRAD * NDOFS(domain)) )
            end if
         end do

      end subroutine ConstructMPIFacesStorage

      subroutine MPI_Face_SendQ(self, domain)
         implicit none
         class(MPI_Face_t)    :: self
         integer, intent(in)  :: domain
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: ierr, dummyreq
   
#ifdef _HAS_MPI_
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_isend(self % Qsend, self % sizeQ, MPI_DOUBLE, domain-1, DEFAULT_TAG, &
                           MPI_COMM_WORLD, dummyreq, ierr)
            call mpi_request_free(dummyreq, ierr)
         end if
#endif

      end subroutine MPI_Face_SendQ

      subroutine MPI_Face_RecvQ(self, domain)
         implicit none
         class(MPI_Face_t)    :: self
         integer, intent(in)  :: domain
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: ierr, dummyreq
   
#ifdef _HAS_MPI_
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_irecv(self % Qrecv, self % sizeQ, MPI_DOUBLE, domain-1, MPI_ANY_TAG, &
                           MPI_COMM_WORLD, self % Qrecv_req, ierr)
         end if
#endif

      end subroutine MPI_Face_RecvQ

      subroutine MPI_Face_SendU_xyz(self, domain)
         implicit none
         class(MPI_Face_t)    :: self
         integer, intent(in)  :: domain
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: ierr, dummyreq
   
#ifdef _HAS_MPI_
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_isend(self % U_xyzsend, self % sizeU_xyz, MPI_DOUBLE, domain-1, &
                           DEFAULT_TAG, MPI_COMM_WORLD, dummyreq, ierr)
            call mpi_request_free(dummyreq, ierr)
         end if
#endif

      end subroutine MPI_Face_SendU_xyz

      subroutine MPI_Face_RecvU_xyz(self, domain)
         implicit none
         class(MPI_Face_t)    :: self
         integer, intent(in)  :: domain
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: ierr, dummyreq
   
#ifdef _HAS_MPI_
         if ( self % no_of_faces .gt. 0 ) then
            call mpi_irecv(self % U_xyzrecv, self % sizeU_xyz, MPI_DOUBLE, domain-1, &
                           DEFAULT_TAG, MPI_COMM_WORLD, self % gradQrecv_req, ierr)
         end if
#endif

      end subroutine MPI_Face_RecvU_xyz

      subroutine MPI_Face_WaitForSolution(self) 
         implicit none 
         class(MPI_Face_t)    :: self
#ifdef _HAS_MPI_ 
! 
!        --------------- 
!        Local variables        
!        --------------- 
! 
         integer              :: ierr 
         integer              :: status(MPI_STATUS_SIZE)
 
!
!        ----------------------------------- 
!        Wait until the solution is received
!        ----------------------------------- 
!
         call mpi_wait(self % Qrecv_req, status, ierr) 
#endif 
         
      end subroutine MPI_Face_WaitForSolution 
 
      subroutine MPI_Face_WaitForGradients(self) 
         implicit none 
         class(MPI_Face_t)    :: self
#ifdef _HAS_MPI_ 
! 
! 
!        --------------- 
!        Local variables        
!        --------------- 
! 
         integer              :: ierr, status(MPI_STATUS_SIZE)
 
!
!        --------------------------------- 
!        Wait until gradients are received
!        -------------------------------- 
!
         call mpi_wait(self % gradQrecv_req, status, ierr) 
#endif 
 
      end subroutine MPI_Face_WaitForGradients 

      subroutine DestructMPIFaces()
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: domain
   
         if ( MPI_Process % doMPIAction ) then
            do domain = 1, MPI_Process % nProcs
               call mpi_faces(domain) % Destruct()
            end do
            safedeallocate(mpi_faces)
         end if

      end subroutine DestructMPIFaces

      type(MPI_Face_t) function Construct_MPI_Face()
         implicit none  

         call MPI_Face_Destruct(Construct_MPI_Face)

      end function Construct_MPI_Face

      subroutine MPI_Face_Construct(self, no_of_faces)
         implicit none
         class(MPI_Face_t) :: self
         integer, intent(in)  :: no_of_faces

         self % no_of_faces = no_of_faces
         allocate(self % faceIDs(no_of_faces))
         allocate(self % elementSide(no_of_faces))

         self % faceIDs       = -1
         self % elementSide   = -1
#ifdef _HAS_MPI_
         self % Qrecv_req     = MPI_REQUEST_NULL
         self % gradQrecv_req = MPI_REQUEST_NULL
#endif

      end subroutine MPI_Face_Construct

      subroutine MPI_Face_Destruct(self)
         implicit none
         class(MPI_Face_t) :: self

         self % no_of_faces = 0
         safedeallocate(self % faceIDs)
         safedeallocate(self % elementSide)
         safedeallocate(self % Qsend)
         safedeallocate(self % U_xyzsend)
         safedeallocate(self % Qrecv)
         safedeallocate(self % U_xyzrecv)
#ifdef _HAS_MPI_
         self % Qrecv_req     = MPI_REQUEST_NULL
         self % gradQrecv_req = MPI_REQUEST_NULL
#endif

      end subroutine MPI_Face_Destruct

end module MPI_Face_Class
