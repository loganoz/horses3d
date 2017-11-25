!
!//////////////////////////////////////////////////////
!
!   @File:    partitioned_mesh.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Sat Nov 25 10:26:09 2017
!   @Last revision date: Sat Nov 25 14:01:20 2017
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 3b34c7a95eb684f3c89837d445c6430f5449f298
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module PartitionedMeshClass
   use SMConstants

   private
   public  PartitionedMesh_t
   public  SendPartitionsMPI, RecvPartitionMPI

   type PartitionedMesh_t
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
   end type PartitionedMesh_t

#ifdef _HAS_MPI_
   type(PartitionedMesh_t), public :: mpi_partition
   type(PartitionedMesh_t), allocatable, public :: mpi_allPartitions(:)
#endif

   interface PartitionedMesh_t
      module procedure  ConstructPartitionedMesh
   end interface

   contains
      subroutine Initialize_MPI_Partitions()
         use MPI_Process_Info
         implicit none
!
!        Create the set of MPI_Partitions in the root rank
!        -------------------------------------------------      
         if ( MPI_Process % doMPIRootAction ) then
            allocate(mpi_allPartitions(MPI_Process % nProcs))
         end if
!
!        Initialize the own MPI partition
!        --------------------------------
         mpi_partition = PartitionedMesh_t(MPI_Process % rank)
         
      end subroutine Initialize_MPI_Partitions
         
      function ConstructPartitionedMesh(ID)
         implicit none
         integer, intent(in)     :: ID
         type(PartitionedMesh_t) :: ConstructPartitionedMesh

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

#endif
      end subroutine RecvPartitionMPI

      subroutine SendPartitionsMPI()
         implicit none
#ifdef _HAS_MPI_

#endif
      end subroutine SendPartitionsMPI
   
end module PartitionedMeshClass
