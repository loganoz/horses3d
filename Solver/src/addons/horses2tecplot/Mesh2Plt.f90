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
         interface
            character(len=LINE_LENGTH) function getFileName( inputLine )
               use SMConstants
               implicit none
               character(len=*)     :: inputLine
            end function getFileName
         end interface
!  
!        ---------------
!        Local variables   
!        ---------------
!
         integer                    :: eID, i, j, k, fid
         type(MeshCoordinates_t)    :: mesh
         character(len=LINE_LENGTH) :: meshPltName
!
!        Read the mesh from the *.hmesh file
!        -----------------------------------
         call mesh % Read(meshFile)
!
!        Create the tecplot mesh file name
!        ---------------------------------
         meshPltName = trim(getFileName(meshFile)) // ".tec"
!
!        Open the file
!        -------------
         open(newunit=fid, file=trim(meshPltName), action="write", status="unknown")
!
!        Write the title and variables
!        -----------------------------
         write(fid,'(A,A,A)') 'TITLE = "',trim(meshFile),'"'
         write(fid,'(A)') 'VARIABLES = "x","y","z"'
!
!        Write the element coordinates
!        -----------------------------
         do eID = 1, mesh % no_of_elements
            associate( e => mesh % elements(eID) )
            write(fid,'(A,I0,A,I0,A,I0,A)') "ZONE I=",e % N(1)+1,", J=",e % N(2)+1, &
                                               ", K=",e % N(3)+1,", F=POINT"
            do k = 0, e % N(3)   ; do j = 0, e % N(2) ;  do i = 0, e % N(1)
               write(fid,'(ES24.16,1X,ES24.16,1X,ES24.16)') e % x(:,i,j,k)
            end do               ; end do             ;  end do
            end associate
         end do
!
!        Close the file
!        --------------
         close(fid)
   
      end subroutine Mesh2Plt

end module Mesh2PltModule
