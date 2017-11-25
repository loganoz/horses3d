!
!//////////////////////////////////////////////////////
!
!   @File:    partitioned_mesh.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Sat Nov 25 10:26:09 2017
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
module MPI_Partitioned_Mesh_MOD
   use SMConstants

   private
   public   MPI_Partitioned_Mesh

   type MPI_Partitioned_Mesh_t
      integer              :: no_of_nodes
      integer              :: no_of_elements
      integer              :: no_of_mpifaces
      integer, allocatable :: nodeIDs(:)
      integer, allocatable :: elementIDs(:)
      integer, allocatable :: mpi_elements(:)
      integer, allocatable :: mpi_element_faceSide(:)
      integer, allocatable :: mpi_face_rotation(:)
      integer, allocatable :: mpi_face_elementSide(:)
      integer, allocatable :: mpi_face_whichProcess(:)
      contains
         procedure   :: Construct => MPI_Partitioned_Mesh_Construct
   end type MPI_Partitioned_Mesh_t
      
   type(MPI_Partitioned_Mesh_t)  :: MPI_Partitioned_Mesh


   contains

      subroutine MPI_Partitioned_Mesh_Construct(self)
         implicit none
         class(MPI_Partitioned_Mesh_t),   intent(out) :: self

      end subroutine MPI_Partitioned_Mesh_Construct
end module MPI_Partitioned_Mesh_MOD
