module ReadMeshFile
   use SMConstants
   use Read_HDF5Mesh_HOPR
   use Read_SpecMesh
   use HexMeshClass
   implicit none

   private
   public constructMeshFromFile, NumOfElemsFromMeshFile

contains
   subroutine constructMeshFromFile( self, fileName, nodes, Nx, Ny, Nz, MeshInnerCurves , dir2D, success )
      implicit none
      !---------------------------------------------------------------
      class(HexMesh)      :: self
      CHARACTER(LEN=*)    :: fileName
      integer             :: nodes
      INTEGER             :: Nx(:), Ny(:), Nz(:)     !<  Polynomial orders for all the elements
      logical             :: MeshInnerCurves         !<  Describe inner curved surfaces? (only for hdf5)
      integer, intent(in) :: dir2D
      LOGICAL             :: success
      !---------------------------------------------------------------
      character(len=LINE_LENGTH) :: ext
      !---------------------------------------------------------------
      interface
         character(len=LINE_LENGTH) function getFileExtension( inputLine )
            use SMConstants
            implicit none
            character(len=*)     :: inputLine
         end function getFileExtension
      end interface
      
      ext = getFileExtension(trim(filename))
      
      if (trim(ext)=='h5') then
         call ConstructMesh_FromHDF5File_( self, fileName, nodes, Nx, Ny, Nz, MeshInnerCurves , dir2D, success )
      elseif (trim(ext)=='mesh') then
         call ConstructMesh_FromSpecMeshFile_( self, fileName, nodes, Nx, Ny, Nz, dir2D, success )
      else
         ERROR STOP 'Mesh file extension not recognized.'
      end if
      
   end subroutine constructMeshFromFile
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function NumOfElemsFromMeshFile(fileName) result(nelem)
      implicit none
      !---------------------------------------------------------------
      CHARACTER(LEN=*)   :: fileName
      integer            :: nelem
      !---------------------------------------------------------------
      character(len=LINE_LENGTH) :: ext
      !---------------------------------------------------------------
      
      interface
         character(len=LINE_LENGTH) function getFileExtension( inputLine )
            use SMConstants
            implicit none
            character(len=*)     :: inputLine
         end function getFileExtension
      end interface
      
      ext = getFileExtension(trim(filename))
      
      if (trim(ext)=='h5') then
         nelem = NumOfElems_HDF5( fileName )
      elseif (trim(ext)=='mesh') then
         nelem = NumOfElems_SpecMesh( fileName )
      else
         ERROR STOP 'Mesh file extension not recognized.'
      end if
      
   end function NumOfElemsFromMeshFile   
end module ReadMeshFile
