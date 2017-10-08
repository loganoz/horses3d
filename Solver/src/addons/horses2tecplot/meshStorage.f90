module meshStorage
   use SMConstants
   use PhysicsStorage
   implicit none
   
   private
   public MeshCoordinates_t

   type ElementCoordinates_t
      integer     :: N(NDIM)
      real(kind=RP), allocatable    :: x(:,:,:,:)
   end type ElementCoordinates_t

   type MeshCoordinates_t
      integer  :: no_of_elements
      type(ElementCoordinates_t),   allocatable    :: elements(:)
      contains
         procedure   :: Read => MeshCoordinates_Read
   end type MeshCoordinates_t

   contains
      subroutine MeshCoordinates_Read(self,meshName)
         use SolutionFile
         implicit none
         class(MeshCoordinates_t)         :: self
         character(len=*), intent(in)     :: meshName
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                          :: arrayDimensions(4)
         integer                          :: fid, eID
         integer  :: i,j,k
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
            e % N(1:3) = arrayDimensions(2:4) - 1 
            allocate( e % x(NDIM,0:e % N(1),0:e % N(2),0:e % N(3)) )
!
!           Read data
!           ---------
            read(fid) e % x

            do k = 0, e % N(3) ; do j = 0, e % N(2);  do i = 0, e % N(1)
               print*, e % x(:,i,j,k)
            end do; end do; end do

            end associate
         end do
!
!        Close file
!        ----------
         close(fid)

      end subroutine MeshCoordinates_Read

      



end module meshStorage
