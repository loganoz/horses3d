!
!//////////////////////////////////////////////////////
!
!   @File:    MeshConsistencySetup.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Dec 24 15:18:47 2017
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
module MeshConsistencySetup
   implicit none
   
   private
   public no_of_meshFiles, meshFileNames

   integer, parameter            :: no_of_meshFiles = 1
   character(len=512), parameter :: meshFileNames(1) = ["../../TestMeshes/CubeWithSphere5.mesh"]

end module MeshConsistencySetup