!
!//////////////////////////////////////////////////////
!
!   @File:    partitioned_mesh.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Sat Nov 25 10:26:09 2017
!   @Last revision date: Mon Nov 27 00:25:05 2017
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 5a757e9f89658144f9cdf9a57dae9d700ca3115a
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module PartitionedMeshClass
   use SMConstants
   use MPI_Process_Info
#ifdef _HAS_MPI_
   use mpi
#endif

   private
   public  PartitionedMesh_t
   
   public  Initialize_MPI_Partitions
   public  SendPartitionsMPI, RecvPartitionMPI

   type PartitionedMesh_t
      logical              :: Constructed
      integer              :: ID
      integer              :: no_of_nodes
      integer              :: no_of_elements
      integer              :: no_of_bdryfaces
      integer, allocatable :: nodeIDs(:)
      integer, allocatable :: elementIDs(:)
      integer, allocatable :: bdryface_elements(:)
      integer, allocatable :: element_bdryfaceSide(:)
      integer, allocatable :: bdryface_rotation(:)
      integer, allocatable :: bdryface_elementSide(:)
      integer, allocatable :: bdryface_sharedDomain(:)
      contains
         procedure   :: Destruct => PartitionedMesh_Destruct
   end type PartitionedMesh_t

   type(PartitionedMesh_t), public :: mpi_partition
   type(PartitionedMesh_t), allocatable, public :: mpi_allPartitions(:)

#ifdef _HAS_MPI_
   integer :: recv_req(8)
   integer, allocatable    :: send_req(:,:)
#endif

   interface PartitionedMesh_t
      module procedure  ConstructPartitionedMesh
   end interface

   contains
      subroutine Initialize_MPI_Partitions()
         implicit none
         integer  :: domain
!
!        Create the set of MPI_Partitions in the root rank
!        -------------------------------------------------      
         if ( MPI_Process % doMPIRootAction ) then
#ifdef _HAS_MPI_
            allocate(mpi_allPartitions(MPI_Process % nProcs))
            allocate(send_req(MPI_Process % nProcs-1,8))
#endif

            do domain = 1, MPI_Process % nProcs
               mpi_allPartitions(domain) = PartitionedMesh_t(domain)
            end do
         end if
!
!        Initialize the own MPI partition
!        --------------------------------
         if ( MPI_Process % doMPIAction ) then
            mpi_partition = PartitionedMesh_t(MPI_Process % rank)
         end if

      end subroutine Initialize_MPI_Partitions
         
      function ConstructPartitionedMesh(ID)
!
!        ********************************************************
!           This is the PartitionedMesh_t constructor
!        ********************************************************
!
         implicit none
         integer, intent(in)     :: ID
         type(PartitionedMesh_t) :: ConstructPartitionedMesh

         ConstructPartitionedMesh % Constructed = .false.
         ConstructPartitionedMesh % ID = ID
         ConstructPartitionedMesh % no_of_nodes = 0
         ConstructPartitionedMesh % no_of_elements = 0
         ConstructPartitionedMesh % no_of_bdryfaces = 0

         safedeallocate(ConstructPartitionedMesh % nodeIDs)
         safedeallocate(ConstructPartitionedMesh % elementIDs)
         safedeallocate(ConstructPartitionedMesh % bdryface_elements)
         safedeallocate(ConstructPartitionedMesh % element_bdryfaceSide)
         safedeallocate(ConstructPartitionedMesh % bdryface_rotation)
         safedeallocate(ConstructPartitionedMesh % bdryface_elementSide)
         safedeallocate(ConstructPartitionedMesh % bdryface_sharedDomain)
   
      end function ConstructPartitionedMesh

      subroutine RecvPartitionMPI()
         implicit none
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: sizes(3), ierr

         if ( MPI_Process % isRoot ) return
!
!        First receive number of nodes, elements, and bdryfaces
!        ------------------------------------------------------
         call mpi_irecv(sizes, 3, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(1), ierr)
!
!        Wait until the message is received
!        ----------------------------------
         call mpi_wait(recv_req(1), MPI_STATUS_IGNORE, ierr)
!
!        Get sizes and allocate
!        ----------------------
         mpi_partition % no_of_nodes = sizes(1)
         mpi_partition % no_of_elements = sizes(2)
         mpi_partition % no_of_bdryfaces = sizes(3)

         allocate(mpi_partition % nodeIDs              (mpi_partition % no_of_nodes    ))
         allocate(mpi_partition % elementIDs           (mpi_partition % no_of_elements ))
         allocate(mpi_partition % bdryface_elements    (mpi_partition % no_of_bdryfaces))
         allocate(mpi_partition % element_bdryfaceSide (mpi_partition % no_of_bdryfaces))
         allocate(mpi_partition % bdryface_rotation    (mpi_partition % no_of_bdryfaces))
         allocate(mpi_partition % bdryface_elementSide (mpi_partition % no_of_bdryfaces))
         allocate(mpi_partition % bdryface_sharedDomain(mpi_partition % no_of_bdryfaces))
!
!        Receive the rest of the PartitionedMesh_t arrays
!        ------------------------------------------------
         call mpi_irecv(mpi_partition % nodeIDs, mpi_partition % no_of_nodes, MPI_INT, 0, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(2), ierr)
   
         call mpi_irecv(mpi_partition % elementIDs, mpi_partition % no_of_elements, MPI_INT, 0, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(3), ierr)

         call mpi_irecv(mpi_partition % bdryface_elements, mpi_partition % no_of_bdryfaces, &
                        MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(4), ierr)

         call mpi_irecv(mpi_partition % element_bdryfaceSide, mpi_partition % no_of_bdryfaces, &
                        MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(5), ierr)

         call mpi_irecv(mpi_partition % bdryface_rotation, mpi_partition % no_of_bdryfaces, &
                        MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(6), ierr)
                     
         call mpi_irecv(mpi_partition % bdryface_elementSide, mpi_partition % no_of_bdryfaces, &
                        MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(7), ierr)

         call mpi_irecv(mpi_partition % bdryface_sharedDomain, mpi_partition % no_of_bdryfaces, &
                        MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(8), ierr)
!
!        Wait until all messages have been received
!        ------------------------------------------
         call mpi_waitall(8, recv_req, MPI_STATUS_IGNORE, ierr) 

         mpi_partition % Constructed = .true.
#endif
      end subroutine RecvPartitionMPI

      subroutine SendPartitionsMPI()
         implicit none
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer          :: sizes(3)
         integer          :: domain, ierr
         integer          :: status(MPI_STATUS_SIZE,MPI_Process % nProcs, 8)
!
!        Send the MPI mesh partition to all processes 
!        --------------------------------------------
         do domain = 2, MPI_Process % nProcs
!
!           Send first the sizes
!           --------------------
            sizes(1) = mpi_allPartitions(domain) % no_of_nodes
            sizes(2) = mpi_allPartitions(domain) % no_of_elements
            sizes(3) = mpi_allPartitions(domain) % no_of_bdryfaces
            call mpi_isend(sizes, 3, MPI_INT, domain-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                           send_req(domain-1,1), ierr)
         end do

         do domain = 2, MPI_Process % nProcs
            call mpi_isend(mpi_allPartitions(domain) % nodeIDs, &
                           mpi_allPartitions(domain) % no_of_nodes, MPI_INT, domain-1, &
                           DEFAULT_TAG, MPI_COMM_WORLD, send_req(domain-1,2), ierr)
   
            call mpi_isend(mpi_allPartitions(domain) % elementIDs, &
                           mpi_allPartitions(domain) % no_of_elements, MPI_INT, domain-1, &
                           DEFAULT_TAG, MPI_COMM_WORLD, send_req(domain-1,3), ierr)

            call mpi_isend(mpi_allPartitions(domain) % bdryface_elements, &
                           mpi_allPartitions(domain) % no_of_bdryfaces, &
                           MPI_INT, domain-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                           send_req(domain-1,4), ierr)

            call mpi_isend(mpi_allPartitions(domain) % element_bdryfaceSide, &
                           mpi_allPartitions(domain) % no_of_bdryfaces, &
                           MPI_INT, domain-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                           send_req(domain-1,5), ierr)

            call mpi_isend(mpi_allPartitions(domain) % bdryface_rotation, &
                           mpi_allPartitions(domain) % no_of_bdryfaces, &
                           MPI_INT, domain-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                           send_req(domain-1,6), ierr)
                     
            call mpi_isend(mpi_allPartitions(domain) % bdryface_elementSide, &
                           mpi_allPartitions(domain) % no_of_bdryfaces, &
                           MPI_INT, domain-1, DEFAULT_TAG, MPI_COMM_WORLD, &  
                           send_req(domain-1,7), ierr)

            call mpi_isend(mpi_allPartitions(domain) % bdryface_sharedDomain, &
                           mpi_allPartitions(domain) % no_of_bdryfaces, &
                           MPI_INT, domain-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                           send_req(domain-1,8), ierr)
         end do
!
!        Copy directly the MPI mesh partition of the root
!        ------------------------------------------------
         mpi_partition = mpi_allPartitions(1)
         mpi_partition % Constructed = .true.
!
!        Wait until all messages have been delivered
!        -------------------------------------------
         call mpi_waitall(8*(MPI_Process % nProcs - 1), send_req, status, ierr) 
!
!        Destruct the partitions
!        -----------------------
         do domain = 1, MPI_Process % nProcs
            call mpi_allPartitions(domain) % Destruct
         end do

#endif
      end subroutine SendPartitionsMPI

      subroutine PartitionedMesh_Destruct(self)
         implicit none
         class(PartitionedMesh_t) :: self

         self % Constructed     = .false.
         self % ID              = 0
         self % no_of_nodes     = 0
         self % no_of_elements  = 0
         self % no_of_bdryfaces = 0

         safedeallocate(self % nodeIDs)
         safedeallocate(self % elementIDs)
         safedeallocate(self % bdryface_elements)
         safedeallocate(self % element_bdryfaceSide)
         safedeallocate(self % bdryface_rotation)
         safedeallocate(self % bdryface_elementSide)
         safedeallocate(self % bdryface_sharedDomain)

      end subroutine PartitionedMesh_Destruct
   
end module PartitionedMeshClass
