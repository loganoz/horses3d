module Mesh2PltModule
   use SMConstants
   use MeshStorage
   implicit none

   private
   public   Mesh2Plt

   contains
      subroutine Mesh2Plt(meshFile)
         implicit none
         character(len=*), intent(in)     :: meshFile
!  
!        ---------------
!        Local variables   
!        ---------------
!
         type(MeshCoordinates_t)     :: mesh

         call mesh % Read(meshFile)
   
      end subroutine Mesh2Plt

end module Mesh2PltModule
