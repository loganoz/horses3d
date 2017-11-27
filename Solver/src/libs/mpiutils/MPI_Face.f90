#include "Includes.h"
module MPI_Face_Class
   implicit none

   private
   public  MPI_Face_t, mpi_faces

   public  ConstructMPIFaces, DestructMPIFaces

   type MPI_Face_t
      integer              :: no_of_faces
      integer, allocatable :: faceIDs(:)
      integer, allocatable :: elementSide(:)
      integer, allocatable :: Qsend_req(:)
      integer, allocatable :: Qrecv_req(:)
      integer, allocatable :: GradQsend_req(:,:)
      integer, allocatable :: GradQrecv_req(:,:)
      contains
         procedure   :: Construct => MPI_Face_Construct
         procedure   :: Destruct => MPI_Face_Destruct
   end type MPI_Face_t

   type(MPI_Face_t), allocatable    :: mpi_faces(:)

   interface MPI_Face_t
      module procedure Construct_MPI_Face
   end interface MPI_Face_t

   contains

      subroutine ConstructMPIFaces()
         use MPI_Process_Info
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

      end subroutine ConstructMPIFaces

      subroutine DestructMPIFaces()
         use MPI_Process_Info
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
         allocate(self % Qsend_req(no_of_faces))
         allocate(self % Qrecv_req(no_of_faces))
         allocate(self % gradQsend_req(3,no_of_faces))
         allocate(self % gradQrecv_req(3,no_of_faces))

         self % faceIDs = -1
         self % elementSide = -1

      end subroutine MPI_Face_Construct

      subroutine MPI_Face_Destruct(self)
         implicit none
         class(MPI_Face_t) :: self

         self % no_of_faces = 0
         safedeallocate(self % faceIDs)
         safedeallocate(self % elementSide)
         safedeallocate(self % Qsend_req)
         safedeallocate(self % Qrecv_req)
         safedeallocate(self % gradQsend_req)
         safedeallocate(self % gradQrecv_req)

      end subroutine MPI_Face_Destruct

end module MPI_Face_Class
