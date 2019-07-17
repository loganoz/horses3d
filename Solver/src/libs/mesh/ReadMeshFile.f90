!
!//////////////////////////////////////////////////////
!
!   @File:    ReadMeshFile.f90
!   @Author:  Andrés Rueda (am.rueda@upm.es)
!   @Created: Sun Apr 27 12:57:00 2017
!   @Last revision date: Wed Jul 17 11:52:46 2019
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 67e046253a62f0e80d1892308486ec5aa1160e53
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module ReadMeshFile
   use SMConstants
   use Read_HDF5Mesh_HOPR
   use Read_SpecMesh
   use HexMeshClass
   use MeshTypes           , only: SPECMESH, HOPRMESH
   use FileReadingUtilities, only: getFileExtension
   implicit none

   private
   public constructMeshFromFile, NumOfElemsFromMeshFile, MeshFileType

contains
   subroutine constructMeshFromFile( self, fileName, nodes, Nx, Ny, Nz, MeshInnerCurves , dir2D, success )
      implicit none
      !---------------------------------------------------------------
      type(HexMesh)                       :: self
      CHARACTER(LEN=*)                    :: fileName
      integer                             :: nodes
      INTEGER                             :: Nx(:), Ny(:), Nz(:)     !<  Polynomial orders for all the elements
      logical                             :: MeshInnerCurves         !<  Describe inner curved surfaces? (only for hdf5)
      integer                             :: dir2D
      LOGICAL           , intent(out)     :: success
      !---------------------------------------------------------------
      character(len=LINE_LENGTH) :: ext
      integer                    :: nelem
      integer                    :: eID
      !---------------------------------------------------------------
      
      ext = getFileExtension(trim(filename))
      
      if (trim(ext)=='h5') then
         call ConstructMesh_FromHDF5File_( self, fileName, nodes, Nx, Ny, Nz, MeshInnerCurves , dir2D, success )
      elseif (trim(ext)=='mesh') then
         call ConstructMesh_FromSpecMeshFile_( self, fileName, nodes, Nx, Ny, Nz, dir2D, success )
      else
         ERROR STOP 'Mesh file extension not recognized.'
      end if
      
      self % NDOF = 0
      do eID=1, self % no_of_elements
         self % NDOF = self % NDOF + product(self % elements(eID) % Nxyz + 1)
      end do
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
      
      ext = getFileExtension(trim(filename))
      
      if (trim(ext)=='h5') then
         nelem = NumOfElems_HDF5( fileName )
      elseif (trim(ext)=='mesh') then
         nelem = NumOfElems_SpecMesh( fileName )
      else
         ERROR STOP 'Mesh file extension not recognized.'
      end if
      
   end function NumOfElemsFromMeshFile  
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   integer function MeshFileType(fileName)
      implicit none
      !-arguments----------------------------------------------------
      CHARACTER(LEN=*)   :: fileName
      !-local-variables----------------------------------------------
      character(len=LINE_LENGTH) :: ext
      !--------------------------------------------------------------
      ext = getFileExtension(trim(filename))
      
      if (trim(ext)=='h5') then
         MeshFileType = HOPRMESH
      elseif (trim(ext)=='mesh') then
         MeshFileType = SPECMESH
      else
         ERROR STOP 'Mesh file extension not recognized.'
      end if
      
   end function MeshFileType
end module ReadMeshFile
