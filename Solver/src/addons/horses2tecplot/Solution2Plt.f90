module Solution2PltModule
   use SMConstants
   use SolutionFile
   implicit none

   private
   public   Solution2Plt

   contains
   
      subroutine Solution2Plt(meshName, solutionName, performInterpolation, Npoints)
         implicit none  
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         logical,          intent(in)     :: performInterpolation
         integer,          intent(in)     :: Npoints
!
!        ---------------
!        Local variables
!        ---------------
!

         if ( performInterpolation ) then

         else
            call Solution2Plt_GaussPoints(meshName, solutionName)

         end if

      end subroutine Solution2Plt

      subroutine Solution2Plt_GaussPoints(meshName, solutionName)
         use MeshStorage
         use SolutionStorage
         implicit none  
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         interface
            character(len=LINE_LENGTH) function getFileName(inputLine)
               use SMConstants
               implicit none
               character(len=*), intent(in)     :: inputLine
            end function getFileName
         end interface
!
!        ---------------
!        Local variables
!        ---------------
!
         type(MeshCoordinates_t)    :: mesh
         type(SolutionStorage_t)    :: solution
         character(len=LINE_LENGTH) :: meshPltName
         character(len=LINE_LENGTH) :: solutionFile
         character(len=1024)        :: title
         integer                    :: no_of_elements, eID
         integer                    :: fid
         integer                    :: Nmesh(4), Nsol(4)

!
!        Read the mesh and solution data
!        -------------------------------
         call mesh % Read(meshName)
         call solution % Read(SolutionName)
!
!        Check that the number of elements in both mesh and solutions are the same
!        -------------------------------------------------------------------------
         if ( mesh % no_of_elements .ne. solution % no_of_elements ) then
            write(STD_OUT,'(A,I0,A,I0,A)') "The number of elements in the mesh (",mesh % no_of_elements,&
                                           ") differs to that of the solution (",solution % no_of_elements,")."
            return
         else
            no_of_elements = mesh % no_of_elements
         end if
!
!        Write the solution file name
!        ----------------------------
         solutionFile = trim(getFileName(solutionName)) // ".tec"
!
!        Create the file
!        ---------------
         open(newunit = fid, file = trim(solutionFile), action = "write", status = "unknown")
!
!        Add the title
!        -------------
         write(title,'(A,A,A,A,A)') '"Generated from ',trim(meshName),' and ',trim(solutionName),'"'
         write(fid,'(A,A)') "TITLE = ", trim(title)
!
!        Add the variables
!        -----------------
         write(fid,'(A)') 'VARIABLES = "x","y","z","rho","rhou","rhov","rhow","rhoe"'
!
!        Write each element zone
!        -----------------------
         do eID = 1, no_of_elements
            associate ( eM => mesh % elements(eID), eS => solution % elements(eID) )
            call WriteElementSolutionGaussPoints(fid,eM % N,eM % x, eS % N, eS % Q)            
            end associate
         end do
!
!        Close the file
!        --------------
         close(fid)
      
      end subroutine Solution2Plt_GaussPoints
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
!     Gauss Points procedures
!     -----------------------
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine WriteElementSolutionGaussPoints(fid, Nmesh, x, Nsol, Q)
         implicit none
         integer,       intent(in) :: fid
         integer,       intent(in) :: Nmesh(3)
         real(kind=RP), intent(in) :: x(1:3,0:Nmesh(1),0:Nmesh(2),0:Nmesh(3))
         integer,       intent(in) :: Nsol(3)
         real(kind=RP), intent(in) :: Q(0:Nsol(1),0:Nsol(2),0:Nsol(3),1:5)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer              :: i,j,k,var
!
!        Write variables: TODO consider p-Adaption
!        ---------------        
         write(fid,'(A,I0,A,I0,A,I0,A)') "ZONE I=",Nmesh(1)+1,", J=",Nmesh(2)+1, &
                                            ", K=",Nmesh(3)+1,", F=POINT"

         do k = 0, Nmesh(3)   ; do j = 0, Nmesh(2)    ; do i = 0, Nmesh(1)
            write(fid,'(ES24.16,1X,ES24.16,1X,ES24.16)',advance="no") x(1,i,j,k), x(2,i,j,k), x(3,i,j,k)
            do var = 1, 5
               write(fid,'(1X,ES24.16)', advance="no") Q(i,j,k,var)
            end do
            write(fid,*)
         end do               ; end do                ; end do

      end subroutine WriteElementSolutionGaussPoints

end module Solution2PltModule
