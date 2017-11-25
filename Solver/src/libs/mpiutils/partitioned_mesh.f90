!
!//////////////////////////////////////////////////////
!
!   @File:    partitioned_mesh.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Sat Nov 25 10:26:09 2017
!   @Last revision date: Sat Nov 25 13:29:58 2017
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 4040089196cab3f09c3e2ea6cb34714311470831
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module PartitionedMeshClass
   use SMConstants

   private
   public  PartitionedMesh_t

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
#endif

   interface PartitionedMesh_t
      module procedure  ConstructPartitionedMesh
   end interface

   contains
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
   
end module PartitionedMeshClass
