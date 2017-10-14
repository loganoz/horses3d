#include "Includes.h"
module Storage
   use SMConstants
   use PhysicsStorage
   implicit none
   
   private
   public Mesh_t, Element_t

   type Element_t
!                                /* Mesh quantities */
      integer                    :: Nmesh(NDIM)
      real(kind=RP), pointer     :: x(:,:,:,:)
!                                /* Solution quantities */
      integer                    :: Nsol(NDIM)
      real(kind=RP), pointer     :: Q(:,:,:,:)
!                                /* Output quantities */
      integer                    :: Nout(NDIM)
      real(kind=RP), pointer     :: xOut(:,:,:,:)
      real(kind=RP), pointer     :: Qout(:,:,:,:)
   end type Element_t

   type Mesh_t
      integer  :: no_of_elements
      type(Element_t),   allocatable    :: elements(:)
      character(len=LINE_LENGTH) :: meshName
      character(len=LINE_LENGTH) :: solutionName
      contains
         procedure   :: ReadMesh     => Mesh_ReadMesh
         procedure   :: ReadSolution => Mesh_ReadSolution
   end type Mesh_t

   contains
      subroutine Mesh_ReadMesh(self,meshName)
         use SolutionFile
         implicit none
         class(Mesh_t)         :: self
         character(len=*), intent(in)     :: meshName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                          :: arrayDimensions(4)
         integer                          :: fid, eID
         integer  :: i,j,k

         self % meshName = trim(meshName)
!
!        Get number of elements
!        ----------------------
         self % no_of_elements = getSolutionFileNoOfElements(meshName)
!
!        Allocate elements
!        -----------------
         allocate(self % elements(self % no_of_elements))
!
!        Read coordinates
!        ----------------
         fid = putSolutionFileInReadDataMode(meshName)
      
         do eID = 1, self % no_of_elements
            associate ( e => self % elements(eID) ) 
            call getSolutionFileArrayDimensions(fid,arrayDimensions)
!
!           Allocate memory for the coordinates
!           -----------------------------------            
            e % Nmesh(1:3) = arrayDimensions(2:4) - 1 
            allocate( e % x(NDIM,0:e % Nmesh(1),0:e % Nmesh(2),0:e % Nmesh(3)) )
!
!           Read data
!           ---------
            read(fid) e % x

            end associate
         end do
!
!        Close file
!        ----------
         close(fid)

      end subroutine Mesh_ReadMesh

      subroutine Mesh_ReadSolution(self,solutionName)
         use SolutionFile
         implicit none
         class(Mesh_t)         :: self
         character(len=*), intent(in)     :: solutionName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                          :: no_of_elements
         integer                          :: arrayDimensions(4)
         integer                          :: fid, eID
         integer  :: i,j,k

         self % solutionName = trim(solutionName)
!
!        Get number of elements
!        ----------------------
         no_of_elements = getSolutionFileNoOfElements(solutionName)
         if ( self % no_of_elements .ne. no_of_elements ) then
            write(STD_OUT,'(A,I0,A,I0,A)') "The number of elements in the mesh (",self % no_of_elements,&
                                           ") differs to that of the solution (",no_of_elements,")."
            errorMessage(STD_OUT)
            stop 
         end if
!
!        Read coordinates
!        ----------------
         fid = putSolutionFileInReadDataMode(solutionName)
      
         do eID = 1, self % no_of_elements
            associate ( e => self % elements(eID) ) 
            call getSolutionFileArrayDimensions(fid,arrayDimensions)
!
!           Allocate memory for the coordinates
!           -----------------------------------            
            e % Nsol(1:3) = arrayDimensions(1:3) - 1 
            allocate( e % Q(0:e % Nsol(1),0:e % Nsol(2),0:e % Nsol(3),1:5) )
!
!           Read data
!           ---------
            read(fid) e % Q

            end associate
         end do
!
!        Close file
!        ----------
         close(fid)

      end subroutine Mesh_ReadSolution

end module Storage
