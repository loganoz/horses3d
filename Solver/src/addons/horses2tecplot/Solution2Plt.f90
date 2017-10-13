module Solution2PltModule
   use SMConstants
   use SolutionFile
   implicit none

   private
   public   Solution2Plt

   type  InterpolationMatrices_t
      logical                          :: Constructed = .false.
      real(kind=RP),    allocatable    :: T(:,:)
   end type InterpolationMatrices_t

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
         integer, parameter      :: Nout(3) = (/10,10,10/)

         if ( performInterpolation ) then

         else
!            call Solution2Plt_GaussPoints(meshName, solutionName)
            call Solution2Plt_GaussPoints_FixedOrder(meshName, solutionName, Nout)

         end if

      end subroutine Solution2Plt
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
!     Gauss Points procedures
!     -----------------------
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
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

      subroutine WriteElementSolutionGaussPoints(fid, Nmesh, x, Nsol, Q, spAM, spAS)
         use NodalStorageClass
         use ProlongMeshAndSolution
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
            call prolongMeshToGaussPoints(Nmesh,x,Nsol,xSol, spAM, spAS)

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
!
!//////////////////////////////////////////////////////////////////////////////////
!
!     Gauss points with fixed order procedures
!     ----------------------------------------
!
!//////////////////////////////////////////////////////////////////////////////////
!
      subroutine Solution2Plt_GaussPoints_FixedOrder(meshName, solutionName, Nout)
         use MeshStorage
         use SolutionStorage
         use NodalStorageClass
         use PolynomialInterpAndDerivsModule
         implicit none  
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         integer,          intent(in)     :: Nout(3)
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
         type(MeshCoordinates_t)                    :: mesh
         type(SolutionStorage_t)                    :: solution
         character(len=LINE_LENGTH)                 :: meshPltName
         character(len=LINE_LENGTH)                 :: solutionFile
         character(len=1024)                        :: title
         integer                                    :: no_of_elements, eID
         integer                                    :: fid
         integer                                    :: Nmesh(4), Nsol(4)
         type(NodalStorage), allocatable            :: spA(:,:,:)
         type(InterpolationMatrices_t), allocatable :: Tset(:,:)
!
!        Allocate nodal storage and interpolation matrices
!        -------------------------------------------------
         allocate( spA(0:NMAX, 0:NMAX, 0:NMAX) )
         allocate( Tset(0:NMAX, 0:NMAX) )
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
!        Allocate the output spectral basis
!        ----------------------------------
         call spA(Nout(1), Nout(2), Nout(3)) % Construct(Nout(1), Nout(2), Nout(3))
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
!           Construct interpolation matrices
!           --------------------------------
            associate( spAS => spA(eS % N(1), eS % N(2), eS % N(3)) , &
                     spAout => spA(Nout(1), Nout(2), Nout(3))           )
            if ( .not. Tset(Nout(1),eS % N(1)) % Constructed ) then
               allocate ( Tset(Nout(1),eS % N(1)) % T(0:Nout(1),0:eS % N(1)) )
               call PolynomialInterpolationMatrix( eS % N(1), Nout(1), spAS % xi, spAS % wbx, spAout % xi, &
                                                   Tset(Nout(1),eS % N(1)) % T)
               Tset(Nout(1),eS % N(1)) % Constructed = .true.
            end if

            if ( .not. Tset(Nout(2),eS % N(2)) % Constructed ) then
               allocate ( Tset(Nout(2),eS % N(2)) % T(0:Nout(2),0:eS % N(2)) )
               call PolynomialInterpolationMatrix( eS % N(2), Nout(2), spAS % xi, spAS % wbx, spAout % xi, &
                                                   Tset(Nout(2),eS % N(2)) % T)
               Tset(Nout(2),eS % N(2)) % Constructed = .true.
            end if

            if ( .not. Tset(Nout(3),eS % N(3)) % Constructed ) then
               allocate ( Tset(Nout(3),eS % N(3)) % T(0:Nout(3),0:eS % N(3)) )
               call PolynomialInterpolationMatrix( eS % N(3), Nout(3), spAS % xi, spAS % wbx, spAout % xi, &
                                                   Tset(Nout(3),eS % N(3)) % T)
               Tset(Nout(3),eS % N(3)) % Constructed = .true.
            end if
            end associate
!
!           Write the tecplot file
!           ----------------------
            call WriteElementSolutionGaussPoints_FixedOrder(fid,eM % N,eM % x, eS % N, eS % Q, Nout, &
                                                                 spA(eM % N(1),eM % N(2),eM % N(3)), & 
                                                                 spA(eS % N(1),eS % N(2),eS % N(3)), &
                                                                       spA(Nout(1),Nout(2),Nout(3)), &
                                                                        Tset(Nout(1),eS % N(1)) % T, &
                                                                        Tset(Nout(2),eS % N(2)) % T, &
                                                                        Tset(Nout(3),eS % N(3)) % T   )
            end associate
         end do
!
!        Close the file
!        --------------
         close(fid)

      end subroutine Solution2Plt_GaussPoints_FixedOrder

      subroutine WriteElementSolutionGaussPoints_FixedOrder(fid, Nmesh, x, Nsol, Q, Nout, spAM, spAS, spAout, Tx, Ty, Tz)
         use NodalStorageClass
         use prolongMeshAndSolution
         implicit none
         integer,            intent(in) :: fid
         integer,            intent(in) :: Nmesh(3)
         real(kind=RP),      intent(in) :: x(1:3,0:Nmesh(1),0:Nmesh(2),0:Nmesh(3))
         integer,            intent(in) :: Nsol(3)
         real(kind=RP),      intent(in) :: Q(0:Nsol(1),0:Nsol(2),0:Nsol(3),1:5)
         integer,            intent(in) :: Nout(3)
         type(NodalStorage), intent(in) :: spAM
         type(NodalStorage), intent(in) :: spAS
         type(NodalStorage), intent(in) :: spAout
         real(kind=RP),      intent(in) :: Tx(0:Nout(1),0:Nsol(1))
         real(kind=RP),      intent(in) :: Ty(0:Nout(2),0:Nsol(2))
         real(kind=RP),      intent(in) :: Tz(0:Nout(3),0:Nsol(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer              :: i,j,k,var
         real(kind=RP)        :: xOut(1:3,0:Nout(1),0:Nout(2),0:Nout(3))
         real(kind=RP)        :: Qout(0:Nout(1),0:Nout(2),0:Nout(3),1:5)
!
!        Project the mesh onto the output polynomial order
!        ---------------------------------------------------
         if ( all( Nmesh .eq. Nout ) ) then
            xOut = x

         else
            call prolongMeshToGaussPoints(Nmesh,x,Nout,xOut, spAM, spAout)

         end if
!
!        Project the solution onto the output polynomial order
!        -----------------------------------------------------
         if ( all( Nsol .eq. Nout ) ) then
            Qout = Q
   
         else
            call prolongSolutionToGaussPoints(Nsol, Q, Nout, Qout, Tx, Ty, Tz)

         end if
!
!        Write variables
!        ---------------        
         write(fid,'(A,I0,A,I0,A,I0,A)') "ZONE I=",Nout(1)+1,", J=",Nout(2)+1, &
                                            ", K=",Nout(3)+1,", F=POINT"

         do k = 0, Nout(3)   ; do j = 0, Nout(2)    ; do i = 0, Nout(1)
            write(fid,'(ES24.16,1X,ES24.16,1X,ES24.16)',advance="no") xOut(1,i,j,k), xOut(2,i,j,k), xOut(3,i,j,k)
            do var = 1, 5
               write(fid,'(1X,ES24.16)', advance="no") Qout(i,j,k,var)
            end do
            write(fid,*)
         end do               ; end do                ; end do

      end subroutine WriteElementSolutionGaussPoints_FixedOrder

end module Solution2PltModule
