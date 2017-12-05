#include "Includes.h"
module MPI_Face_Class
   implicit none

   private
   public  MPI_Face_t, mpi_faces, MPI_Faces_Constructed

   public  ConstructMPIFaces, DestructMPIFaces
   public  WaitUntilSolutionIsReady 
   public  WaitUntilGradientsAreReady 

   type MPI_Face_t
      integer              :: no_of_faces
      integer, allocatable :: faceIDs(:)
      integer, allocatable :: elementSide(:)
      integer, allocatable :: Qrecv_req(:)
      integer, allocatable :: GradQrecv_req(:,:)
      contains
         procedure   :: Construct => MPI_Face_Construct
         procedure   :: Destruct => MPI_Face_Destruct
   end type MPI_Face_t

   logical                          :: MPI_Faces_Constructed = .false.
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

         MPI_Faces_Constructed = .true.

      end subroutine ConstructMPIFaces

      subroutine WaitUntilSolutionIsReady(faces) 
         use MPI_Process_Info 
#ifdef _HAS_MPI_ 
         use mpi 
#endif 
         implicit none 
         type(MPI_Face_t)     :: faces(MPI_Process % nProcs) 
! 
!        --------------- 
!        Local variables        
!        --------------- 
! 
         integer              :: domain, ierr 
         integer              :: no_of_faces, lb, ub
         integer, allocatable :: status(:,:), requests(:) 
 
#ifdef _HAS_MPI_ 
!
!        ---------------------------------
!        Get the total number of MPI Faces 
!        ---------------------------------
!
         no_of_faces = 0
         do domain = 1, MPI_Process % nProcs 
            no_of_faces = no_of_faces + faces(domain) % no_of_faces
         end do

         allocate(status(MPI_STATUS_SIZE, no_of_faces))
         allocate(requests(no_of_faces))
!
!        ------------------------------
!        Get the requests in a 1D array
!        ------------------------------
!
         no_of_faces = 0
         do domain = 1, MPI_Process % nProcs
            if ( faces(domain) % no_of_faces .le. 0 ) cycle
!
!           Lower and upper bounds
!           ----------------------
            lb = no_of_faces + 1
            ub = lb - 1 + faces(domain) % no_of_faces
!
!           Gather requests
!           ---------------
            requests(lb:ub) = faces(domain) % Qrecv_req
!
!           Get next face
!           -------------
            no_of_faces = no_of_faces + faces(domain) % no_of_faces

         end do
!
!        ----------------------------------- 
!        Wait until the solution is received
!        ----------------------------------- 
!
         do domain = 1, MPI_Process % nProcs
            call mpi_waitall(no_of_faces, requests, status, ierr) 
         end do
#endif 
         
      end subroutine WaitUntilSolutionIsReady 
 
      subroutine WaitUntilGradientsAreReady(faces) 
         use MPI_Process_Info 
#ifdef _HAS_MPI_ 
         use mpi 
#endif 
         implicit none 
         type(MPI_Face_t)     :: faces(MPI_Process % nProcs) 
! 
! 
!        --------------- 
!        Local variables        
!        --------------- 
! 
         integer              :: domain, ierr, i, j, lb, ub
         integer              :: no_of_faces
         integer, allocatable :: status(:,:), requests(:) 
 
#ifdef _HAS_MPI_ 
!
!        ---------------------------------
!        Get the total number of MPI Faces 
!        ---------------------------------
!
         no_of_faces = 0
         do domain = 1, MPI_Process % nProcs 
            no_of_faces = no_of_faces + faces(domain) % no_of_faces
         end do

         allocate(status(MPI_STATUS_SIZE, 3*no_of_faces))
         allocate(requests(3*no_of_faces))
!
!        ------------------------------
!        Get the requests in a 1D array
!        ------------------------------
!
         no_of_faces = 0
         do domain = 1, MPI_Process % nProcs
            if ( faces(domain) % no_of_faces .eq. 0 ) cycle
            do j = 1, 3
!
!           Lower and upper bounds
!           ----------------------
            lb = 3*no_of_faces + 1 + (j-1)*faces(domain) % no_of_faces
            ub = lb - 1 + faces(domain) % no_of_faces
!
!           Gather requests
!           ---------------
            requests(lb:ub) = faces(domain) % gradQrecv_req(j,:)

            end do
!
!           Get next face
!           -------------
            no_of_faces = no_of_faces + faces(domain) % no_of_faces
         end do
!
!        ----------------------------------- 
!        Wait until the solution is received
!        ----------------------------------- 
!
         call mpi_waitall(3*no_of_faces, requests, status, ierr) 

         deallocate(status) 
         deallocate(requests)
#endif 
 
      end subroutine WaitUntilGradientsAreReady 

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
         allocate(self % Qrecv_req(no_of_faces))
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
         safedeallocate(self % Qrecv_req)
         safedeallocate(self % gradQrecv_req)

      end subroutine MPI_Face_Destruct

end module MPI_Face_Class
