!
!//////////////////////////////////////////////////////
!
!   @File:    MeshPartitioning.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Sat Nov 25 10:26:08 2017
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
module MeshPartitioning
   use SMConstants
   use HexMeshClass

   private
   public   PerformMeshPartitioning


   contains
      subroutine PerformMeshPartitioning(mesh, no_of_domains)
         type(HexMesh), intent(in)  :: mesh
         integer,       intent(in)  :: no_of_domains
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: elementDomain(mesh % no_of_elements)
!
!        ************************************************
!        Now, they will just be ordered in a hard-coded
!        way, and just valid for the TaylorGreen geometry
!        It is required to consider using METIS or a 
!        universal partitioner later.
!        ************************************************
!




      end subroutine PerformMeshPartitioning

      subroutine GetElementsDomain(mesh)
         implicit none
         type(HexMesh), intent(in)  :: mesh

      end subroutine GetElementsDomain

end module MeshPartitioning
