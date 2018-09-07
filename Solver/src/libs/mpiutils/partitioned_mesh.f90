!
!//////////////////////////////////////////////////////
!
!   @File:    partitioned_mesh.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Sat Nov 25 10:26:09 2017
!   @Last revision date: Thu Sep  6 15:28:01 2018
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: bbb1eccef1b477a3cf37da8b42ada5792b1c1bf3
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
      integer              :: no_of_mpifaces
      integer, allocatable :: nodeIDs(:)
      integer, allocatable :: HOPRnodeIDs(:)
      integer, allocatable :: elementIDs(:)
      integer, allocatable :: mpiface_elements(:)
      integer, allocatable :: element_mpifaceSide(:)
      integer, allocatable :: mpiface_rotation(:)
      integer, allocatable :: mpiface_elementSide(:)
      integer, allocatable :: mpiface_sharedDomain(:)
      contains
         procedure   :: Destruct => PartitionedMesh_Destruct
   end type PartitionedMesh_t

   type(PartitionedMesh_t), public :: mpi_partition
   type(PartitionedMesh_t), allocatable, public :: mpi_allPartitions(:)
   
   integer, protected, public :: MPI_Partitioning
   integer, parameter, public :: METIS_PARTITIONING = 1
   integer, parameter, public :: SFC_PARTITIONING   = 2

#ifdef _HAS_MPI_
   integer :: recv_req(8)
   integer, allocatable    :: send_req(:,:)
#endif

   interface PartitionedMesh_t
      module procedure  ConstructPartitionedMesh
   end interface

   contains
      subroutine Initialize_MPI_Partitions(partitioning)
         implicit none
         !-arguments----------------------------------------------------
         character(len=*), intent(in)  :: partitioning
         !-local-variables----------------------------------------------
         integer  :: domain
         !--------------------------------------------------------------
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
            
!           Partitioning method
!           -------------------
            select case (partitioning)
!     
!              Space-filling curve partitioning
!              --------------------------------
               case ('SFC')
                  MPI_Partitioning = SFC_PARTITIONING
!     
!              METIS partitioning
!              ------------------
               case default
                  MPI_Partitioning = METIS_PARTITIONING
#ifndef _HAS_METIS_
                  write(STD_OUT,*) 'ERROR: Metis not linked properly. Two options available:'
                  write(STD_OUT,*) '   * Set the METIS_HOME variable and recompile'
                  write(STD_OUT,*) '   * Use "partitioning = SFC" if the mesh elements are numbered as a space-filling curve'
                  stop
#endif
            end select
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
         ConstructPartitionedMesh % no_of_mpifaces = 0

         safedeallocate(ConstructPartitionedMesh % nodeIDs)
         safedeallocate(ConstructPartitionedMesh % elementIDs)
         safedeallocate(ConstructPartitionedMesh % mpiface_elements)
         safedeallocate(ConstructPartitionedMesh % element_mpifaceSide)
         safedeallocate(ConstructPartitionedMesh % mpiface_rotation)
         safedeallocate(ConstructPartitionedMesh % mpiface_elementSide)
         safedeallocate(ConstructPartitionedMesh % mpiface_sharedDomain)
   
      end function ConstructPartitionedMesh

      subroutine RecvPartitionMPI(meshIsHOPR)
         implicit none
!
!        ---------
!        Arguments
!        ---------
!
         logical, intent(in) :: meshIsHOPR
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: sizes(3), ierr, array_of_statuses(MPI_STATUS_SIZE, 8), status_HOPR(MPI_STATUS_SIZE)
         integer  :: recv_reqHOPR

         if ( MPI_Process % isRoot ) return
!
!        First receive number of nodes, elements, and mpifaces
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
         mpi_partition % no_of_mpifaces = sizes(3)

         allocate(mpi_partition % nodeIDs             (mpi_partition % no_of_nodes   )) 
         allocate(mpi_partition % elementIDs          (mpi_partition % no_of_elements))
         allocate(mpi_partition % mpiface_elements    (mpi_partition % no_of_mpifaces))
         allocate(mpi_partition % element_mpifaceSide (mpi_partition % no_of_mpifaces))
         allocate(mpi_partition % mpiface_rotation    (mpi_partition % no_of_mpifaces))
         allocate(mpi_partition % mpiface_elementSide (mpi_partition % no_of_mpifaces))
         allocate(mpi_partition % mpiface_sharedDomain(mpi_partition % no_of_mpifaces))
         
         if (meshIsHOPR) allocate(mpi_partition % HOPRnodeIDs(mpi_partition % no_of_nodes   )) 
!
!        Receive the rest of the PartitionedMesh_t arrays
!        ------------------------------------------------
         call mpi_irecv(mpi_partition % nodeIDs, mpi_partition % no_of_nodes, MPI_INT, 0, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(2), ierr)
         
         call mpi_irecv(mpi_partition % elementIDs, mpi_partition % no_of_elements, MPI_INT, 0, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(3), ierr)

         call mpi_irecv(mpi_partition % mpiface_elements, mpi_partition % no_of_mpifaces, &
                        MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(4), ierr)

         call mpi_irecv(mpi_partition % element_mpifaceSide, mpi_partition % no_of_mpifaces, &
                        MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(5), ierr)

         call mpi_irecv(mpi_partition % mpiface_rotation, mpi_partition % no_of_mpifaces, &
                        MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(6), ierr)
                     
         call mpi_irecv(mpi_partition % mpiface_elementSide, mpi_partition % no_of_mpifaces, &
                        MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(7), ierr)

         call mpi_irecv(mpi_partition % mpiface_sharedDomain, mpi_partition % no_of_mpifaces, &
                        MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(8), ierr)
         
         if (meshIsHOPR) then
            call mpi_irecv(mpi_partition % HOPRnodeIDs, mpi_partition % no_of_nodes, MPI_INT, 0, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, recv_reqHOPR, ierr)
         end if
!
!        Wait until all messages have been received
!        ------------------------------------------
         call mpi_waitall(8, recv_req, array_of_statuses, ierr) 
         if (meshIsHOPR) call mpi_wait(recv_reqHOPR, status_HOPR, ierr) 

         mpi_partition % Constructed = .true.
#endif
      end subroutine RecvPartitionMPI

      subroutine SendPartitionsMPI(meshIsHOPR)
         implicit none
!
!        ---------
!        Arguments
!        ---------
!
         logical, intent(in) :: meshIsHOPR
#ifdef _HAS_MPI_
!
!        ---------------
!        Local variables
!        ---------------
!
         integer          :: sizes(3)
         integer          :: domain, ierr, msg
         integer          :: array_of_statuses(MPI_STATUS_SIZE,MPI_Process % nProcs)
         integer          :: send_reqHOPR( MPI_Process % nProcs - 1)
!
!        Send the MPI mesh partition to all processes 
!        --------------------------------------------
         do domain = 2, MPI_Process % nProcs
!
!           Send first the sizes
!           --------------------
            sizes(1) = mpi_allPartitions(domain) % no_of_nodes
            sizes(2) = mpi_allPartitions(domain) % no_of_elements
            sizes(3) = mpi_allPartitions(domain) % no_of_mpifaces
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

            call mpi_isend(mpi_allPartitions(domain) % mpiface_elements, &
                           mpi_allPartitions(domain) % no_of_mpifaces, &
                           MPI_INT, domain-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                           send_req(domain-1,4), ierr)

            call mpi_isend(mpi_allPartitions(domain) % element_mpifaceSide, &
                           mpi_allPartitions(domain) % no_of_mpifaces, &
                           MPI_INT, domain-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                           send_req(domain-1,5), ierr)

            call mpi_isend(mpi_allPartitions(domain) % mpiface_rotation, &
                           mpi_allPartitions(domain) % no_of_mpifaces, &
                           MPI_INT, domain-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                           send_req(domain-1,6), ierr)
                     
            call mpi_isend(mpi_allPartitions(domain) % mpiface_elementSide, &
                           mpi_allPartitions(domain) % no_of_mpifaces, &
                           MPI_INT, domain-1, DEFAULT_TAG, MPI_COMM_WORLD, &  
                           send_req(domain-1,7), ierr)

            call mpi_isend(mpi_allPartitions(domain) % mpiface_sharedDomain, &
                           mpi_allPartitions(domain) % no_of_mpifaces, &
                           MPI_INT, domain-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                           send_req(domain-1,8), ierr)
            
            if (meshIsHOPR) then
               call mpi_isend(mpi_allPartitions(domain) % HOPRnodeIDs, &
                           mpi_allPartitions(domain) % no_of_nodes, MPI_INT, domain-1, &
                           DEFAULT_TAG, MPI_COMM_WORLD, send_reqHOPR(domain-1), ierr)
            end if
         end do
!
!        Copy directly the MPI mesh partition of the root
!        ------------------------------------------------
         mpi_partition = mpi_allPartitions(1)
         mpi_partition % Constructed = .true.
!
!        Wait until all messages have been delivered
!        -------------------------------------------
         do msg = 1, 8
            call mpi_waitall(MPI_Process % nProcs - 1, send_req(:,msg), array_of_statuses, ierr) 
         end do
         if (meshIsHOPR) call mpi_waitall(MPI_Process % nProcs - 1, send_reqHOPR(:), array_of_statuses, ierr) 
         
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
         self % no_of_mpifaces = 0

         safedeallocate(self % nodeIDs             )
         safedeallocate(self % HOPRnodeIDs         )
         safedeallocate(self % elementIDs          )
         safedeallocate(self % mpiface_elements    )
         safedeallocate(self % element_mpifaceSide )
         safedeallocate(self % mpiface_rotation    )
         safedeallocate(self % mpiface_elementSide )
         safedeallocate(self % mpiface_sharedDomain)
         

      end subroutine PartitionedMesh_Destruct
   
end module PartitionedMeshClass
