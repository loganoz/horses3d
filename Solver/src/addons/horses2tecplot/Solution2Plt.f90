module Solution2PltModule
   use SMConstants
   use SolutionFile
   implicit none

   private
   public   Solution2Plt

   integer, parameter   :: NMAX = 40

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
         use NodalStorageClass
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
         type(MeshCoordinates_t)         :: mesh
         type(SolutionStorage_t)         :: solution
         character(len=LINE_LENGTH)      :: meshPltName
         character(len=LINE_LENGTH)      :: solutionFile
         character(len=1024)             :: title
         integer                         :: no_of_elements, eID
         integer                         :: fid
         integer                         :: Nmesh(4), Nsol(4)
         type(NodalStorage), allocatable :: spA(:,:,:)

!
!        Read the mesh and solution data
!        -------------------------------
         call mesh % Read(meshName)
         call solution % Read(SolutionName)
!
!        Create the nodal approximation structure
!        ----------------------------------------
         allocate( spA(0:NMAX,0:NMAX,0:NMAX) )
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
!
!           Construct spectral basis
!           ------------------------
            if ( .not. spA(eM % N(1),eM % N(2), eM % N(3)) % Constructed ) then
               call spA(eM % N(1), eM % N(2), eM % N(3) ) % Construct( eM % N(1), eM % N(2), eM % N(3) )
            end if

            if ( .not. spA(eS % N(1),eS % N(2), eS % N(3)) % Constructed ) then
               call spA(eS % N(1), eS % N(2), eS % N(3) ) % Construct( eS % N(1), eS % N(2), eS % N(3) )
            end if
!
!           Write the tecplot file
!           ----------------------
            call WriteElementSolutionGaussPoints(fid,eM % N,eM % x, eS % N, eS % Q, &
                                                spA(eM % N(1),eM % N(2),eM % N(3)), & 
                                                spA(eS % N(1),eS % N(2),eS % N(3))  )
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
      subroutine WriteElementSolutionGaussPoints(fid, Nmesh, x, Nsol, Q, spAM, spAS)
         use NodalStorageClass
         use MeshInterpolation
         implicit none
         integer,            intent(in) :: fid
         integer,            intent(in) :: Nmesh(3)
         real(kind=RP),      intent(in) :: x(1:3,0:Nmesh(1),0:Nmesh(2),0:Nmesh(3))
         integer,            intent(in) :: Nsol(3)
         real(kind=RP),      intent(in) :: Q(0:Nsol(1),0:Nsol(2),0:Nsol(3),1:5)
         type(NodalStorage), intent(in) :: spAM
         type(NodalStorage), intent(in) :: spAS
!
!        ---------------
!        Local variables
!        ---------------
!
         integer              :: i,j,k,var
         real(kind=RP)        :: xSol(1:3,0:Nsol(1),0:Nsol(2),0:NSol(3))
!
!        Project the mesh onto the solution polynomial order
!        ---------------------------------------------------
         if ( all( Nmesh .eq. Nsol ) ) then
            xSol = x

         else
            call prolongMeshToSolutionOrder(Nmesh,x,Nsol,xSol, spAM, spAS)

         end if
!
!        Write variables
!        ---------------        
         write(fid,'(A,I0,A,I0,A,I0,A)') "ZONE I=",Nsol(1)+1,", J=",Nsol(2)+1, &
                                            ", K=",Nsol(3)+1,", F=POINT"

         do k = 0, Nsol(3)   ; do j = 0, Nsol(2)    ; do i = 0, Nsol(1)
            write(fid,'(ES24.16,1X,ES24.16,1X,ES24.16)',advance="no") xSol(1,i,j,k), xSol(2,i,j,k), xSol(3,i,j,k)
            do var = 1, 5
               write(fid,'(1X,ES24.16)', advance="no") Q(i,j,k,var)
            end do
            write(fid,*)
         end do               ; end do                ; end do

      end subroutine WriteElementSolutionGaussPoints

end module Solution2PltModule
