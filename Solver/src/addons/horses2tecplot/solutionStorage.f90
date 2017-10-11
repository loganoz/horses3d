module solutionStorage
   use SMConstants
   use PhysicsStorage
   implicit none
   
   private
   public SolutionStorage_t

   type ElementStorage_t
      integer     :: N(NDIM)
      real(kind=RP), allocatable    :: Q(:,:,:,:)
   end type ElementStorage_t

   type SolutionStorage_t
      integer  :: no_of_elements
      type(ElementStorage_t),   allocatable    :: elements(:)
      contains
         procedure   :: Read => SolutionStorage_Read
   end type SolutionStorage_t

   contains
      subroutine SolutionStorage_Read(self,solutionName)
         use SolutionFile
         implicit none
         class(SolutionStorage_t)         :: self
         character(len=*), intent(in)     :: solutionName
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
         self % no_of_elements = getSolutionFileNoOfElements(solutionName)
!
!        Allocate elements
!        -----------------
         allocate(self % elements(self % no_of_elements))
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
            e % N(1:3) = arrayDimensions(1:3) - 1 
            allocate( e % Q(0:e % N(1),0:e % N(2),0:e % N(3),1:5) )
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

      end subroutine SolutionStorage_Read

end module solutionStorage
