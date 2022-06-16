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
