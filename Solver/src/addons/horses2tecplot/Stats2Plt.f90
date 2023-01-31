#include "Includes.h"
module Stats2PltModule
   use SMConstants
   use SolutionFile
   use Headers
   use InterpolationMatrices 
   use FileReadingUtilities      , only: getFileName
   use Solution2PltModule        , only: WriteBoundaryToTecplot
   implicit none

   private
   public   Stats2Plt

#define PRECISION_FORMAT "(E13.5)"

   contains
      subroutine Stats2Plt(meshName, solutionName, fixedOrder, basis, Nout)
         use getTask
         implicit none  
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         integer,          intent(in)     :: basis
         logical,          intent(in)     :: fixedOrder
         integer,          intent(in)     :: Nout(3)
   
         write(STD_OUT,'(/)')
         call SubSection_Header("Job description")

         select case ( basis )

         case(EXPORT_GAUSS)

            if ( fixedOrder ) then
               write(STD_OUT,'(30X,A3,A)') "->", " Export to Gauss points with fixed order"
               write(STD_OUT,'(30X,A,A30,I0,A,I0,A,I0,A)') "->" , "Output order: [",&
                                                Nout(1),",",Nout(2),",",Nout(3),"]."
               call Stats2Plt_GaussPoints_FixedOrder(meshName, solutionName, Nout)
   
            else
               write(STD_OUT,'(30X,A3,A)') "->", " Export to Gauss points"
               call Stats2Plt_GaussPoints(meshName, solutionName)

            end if

         case(EXPORT_HOMOGENEOUS)
            
            write(STD_OUT,'(30X,A3,A)') "->", " Export to homogeneous points"
            write(STD_OUT,'(30X,A,A30,I0,A,I0,A,I0,A)') "->" , "Output order: [",&
                                        Nout(1),",",Nout(2),",",Nout(3),"]."
            call Stats2Plt_Homogeneous(meshName, solutionName, Nout)

         end select

      end subroutine Stats2Plt
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
!     Gauss Points procedures
!     -----------------------
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Stats2Plt_GaussPoints(meshName, solutionName)
         use Storage
         use NodalStorageClass
         use SharedSpectralBasis
         use OutputVariables
         implicit none  
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
!
!        ---------------
!        Local variables
!        ---------------
!
         type(Mesh_t)                    :: mesh
         character(len=LINE_LENGTH)      :: meshPltName
         character(len=LINE_LENGTH)      :: solutionFile
         character(len=1024)             :: title
         integer                         :: no_of_elements, eID
         integer                         :: fid, bID
         integer                         :: Nmesh(4), Nsol(4)
!
!        Read the mesh and solution data
!        -------------------------------
         call mesh % ReadMesh(meshName)
         call mesh % ReadSolution(SolutionName)
         no_of_elements = mesh % no_of_elements
!
!        Transform zones to the output variables
!        ---------------------------------------
         do eID = 1, no_of_elements
            associate ( e => mesh % elements(eID) )
            e % Nout = e % Nsol
!
!           Construct spectral basis
!           ------------------------
            call addNewSpectralBasis(spA, e % Nmesh, mesh % nodeType)
            call addNewSpectralBasis(spA, e % Nsol, mesh % nodeType)
!
!           Project mesh and solution
!           -------------------------
            call ProjectStorageGaussPoints(e, spA, e % Nmesh, e % Nsol)

            end associate
         end do
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
         write(fid,'(A,A)') 'VARIABLES = "x","y","z","Umean","Vmean","Wmean","Sxx","Syy","Szz","Sxy","Sxz","Syz"'
!
!        Write each element zone
!        -----------------------
         do eID = 1, no_of_elements
            associate ( e => mesh % elements(eID) )
!
!           Write the tecplot file
!           ----------------------
            call WriteElementToTecplot(fid, e, mesh % refs) 
            end associate
         end do
!
!        Write boundaries
!        ----------------
         if (hasBoundaries) then
            do bID=1, size (mesh % boundaries)
               call WriteBoundaryToTecplot(fid, mesh % boundaries(bID), mesh % elements)
            end do
         end if
!
!        Close the file
!        --------------
         close(fid)
      
      end subroutine Stats2Plt_GaussPoints

      subroutine ProjectStorageGaussPoints(e, spA, N1, N2)
         use Storage
         use NodalStorageClass
         use ProlongMeshAndSolution
         implicit none
         type(Element_t)     :: e
         type(NodalStorage_t), intent(in) :: spA(0:)
         integer           , intent(in) :: N1(3)
         integer           , intent(in) :: N2(3)
         
         e % Nout = e % Nsol
         if ( all(e % Nmesh .eq. e % Nout) ) then
            e % xOut(1:,0:,0:,0:) => e % x

         else
            allocate( e % xOut(1:3,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            call prolongMeshToGaussPoints(e, spA, N1, N2)

         end if

         e % statsout(1:,0:,0:,0:) => e % stats

      end subroutine ProjectStorageGaussPoints
!
!//////////////////////////////////////////////////////////////////////////////////
!
!     Gauss points with fixed order procedures
!     ----------------------------------------
!
!//////////////////////////////////////////////////////////////////////////////////
!
      subroutine Stats2Plt_GaussPoints_FixedOrder(meshName, solutionName, Nout)
         use Storage
         use NodalStorageClass
         use SharedSpectralBasis
         use OutputVariables
         implicit none  
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         integer,          intent(in)     :: Nout(3)
!
!        ---------------
!        Local variables
!        ---------------
!
         type(Mesh_t)                               :: mesh
         character(len=LINE_LENGTH)                 :: meshPltName
         character(len=LINE_LENGTH)                 :: solutionFile
         character(len=1024)                        :: title
         integer                                    :: no_of_elements, eID
         integer                                    :: fid, bID
!
!        Read the mesh and solution data
!        -------------------------------
         call mesh % ReadMesh(meshName)
         call mesh % ReadSolution(SolutionName)
!
!        Allocate the output spectral basis
!        ----------------------------------
         call spA(Nout(1)) % Construct(GAUSS, Nout(1))
         call spA(Nout(2)) % Construct(GAUSS, Nout(2))
         call spA(Nout(3)) % Construct(GAUSS, Nout(3))
!
!        Write each element zone
!        -----------------------
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            e % Nout = Nout
!
!           Construct spectral basis
!           ------------------------
            call addNewSpectralBasis(spA, e % Nmesh, mesh % nodeType)
            call addNewSpectralBasis(spA, e % Nsol , mesh % nodeType)
!
!           Construct interpolation matrices
!           --------------------------------
            associate( spAoutXi   => spA(Nout(1)), &
                       spAoutEta  => spA(Nout(2)), &
                       spAoutZeta => spA(Nout(3)) )
            call addNewInterpolationMatrix(Tset, e % Nsol(1), spA(e % Nsol(1)), e % Nout(1), spAoutXi   % x)
            call addNewInterpolationMatrix(Tset, e % Nsol(2), spA(e % Nsol(2)), e % Nout(2), spAoutEta  % x)
            call addNewInterpolationMatrix(Tset, e % Nsol(3), spA(e % Nsol(3)), e % Nout(3), spAoutZeta % x)
            end associate
!
!           Perform interpolation
!           ---------------------
            call ProjectStorageGaussPoints_FixedOrder(e, spA, e % Nmesh, e % Nsol, e % Nout, &
                                                                    Tset(e % Nout(1), e % Nsol(1)) % T, &
                                                                    Tset(e % Nout(2), e % Nsol(2)) % T, &
                                                                    Tset(e % Nout(3), e % Nsol(3)) % T    )

            end associate
         end do
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
         write(fid,'(A,A)') 'VARIABLES = "x","y","z","Umean","Vmean","Wmean","Sxx","Syy","Szz","Sxy","Sxz","Syz"'
!
!        Write elements
!        --------------
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )

            call WriteElementToTecplot(fid, e, mesh % refs)
            end associate
         end do
!
!        Write boundaries
!        ----------------
         if (hasBoundaries) then
            do bID=1, size (mesh % boundaries)
               call WriteBoundaryToTecplot(fid, mesh % boundaries(bID), mesh % elements)
            end do
         end if

!
!        Close the file
!        --------------
         close(fid)

      end subroutine Stats2Plt_GaussPoints_FixedOrder

      subroutine ProjectStorageGaussPoints_FixedOrder(e, spA, NM, NS, Nout, Tx, Ty, Tz)
         use Storage
         use NodalStorageClass
         use ProlongMeshAndSolution
         implicit none
         type(Element_t)     :: e
         type(NodalStorage_t),  intent(in)  :: spA(0:)
         integer           ,  intent(in)  :: NM(3)
         integer           ,  intent(in)  :: NS(3)
         integer           ,  intent(in)  :: Nout(3)
         real(kind=RP),       intent(in)  :: Tx(0:e % Nout(1), 0:e % Nsol(1))
         real(kind=RP),       intent(in)  :: Ty(0:e % Nout(2), 0:e % Nsol(2))
         real(kind=RP),       intent(in)  :: Tz(0:e % Nout(3), 0:e % Nsol(3))
!
!        Project mesh
!        ------------         
         if ( all(e % Nmesh .eq. e % Nout) ) then
            e % xOut(1:,0:,0:,0:) => e % x

         else
            allocate( e % xOut(1:3,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            call prolongMeshToGaussPoints(e, spA, NM, Nout)

         end if
!
!        Project the solution
!        --------------------
         if ( all( e % Nsol .eq. e % Nout ) ) then
            e % statsout(1:,0:,0:,0:) => e % stats
   
         else
            allocate( e % statsout(1:9,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            call prolongSolutionToGaussPoints(9, e % Nsol, e % stats, e % Nout, e % statsout, Tx, Ty, Tz)

         end if

      end subroutine ProjectStorageGaussPoints_FixedOrder
!
!////////////////////////////////////////////////////////////////////////////
!
!     Homogeneous procedures
!     ----------------------
!
!////////////////////////////////////////////////////////////////////////////
!
      subroutine Stats2Plt_Homogeneous(meshName, solutionName, Nout)
         use Storage
         use NodalStorageClass
         use SharedSpectralBasis
         use OutputVariables
         implicit none  
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         integer,          intent(in)     :: Nout(3)
!
!        ---------------
!        Local variables
!        ---------------
!
         type(Mesh_t)                               :: mesh
         character(len=LINE_LENGTH)                 :: meshPltName
         character(len=LINE_LENGTH)                 :: solutionFile
         character(len=1024)                        :: title
         integer                                    :: no_of_elements, eID
         integer                                    :: fid, bID
         real(kind=RP)                              :: xi(0:Nout(1)), eta(0:Nout(2)), zeta(0:Nout(3))
         integer                                    :: i
!
!        Read the mesh and solution data
!        -------------------------------
         call mesh % ReadMesh(meshName)
         call mesh % ReadSolution(SolutionName)
!
!        Set homogeneous nodes
!        ---------------------
         xi   = RESHAPE( (/ (-1.0_RP + 2.0_RP*i/Nout(1),i=0,Nout(1)) /), (/ Nout(1)+1 /) )
         eta  = RESHAPE( (/ (-1.0_RP + 2.0_RP*i/Nout(2),i=0,Nout(2)) /), (/ Nout(2)+1 /) )
         zeta = RESHAPE( (/ (-1.0_RP + 2.0_RP*i/Nout(3),i=0,Nout(3)) /), (/ Nout(3)+1 /) )
!
!        Write each element zone
!        -----------------------
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            e % Nout = Nout
!
!           Construct spectral basis for both mesh and solution
!           ---------------------------------------------------
            call addNewSpectralBasis(spA, e % Nmesh, mesh % nodeType)
            call addNewSpectralBasis(spA, e % Nsol , mesh % nodeType)
!
!           Construct interpolation matrices for the mesh
!           ---------------------------------------------
            call addNewInterpolationMatrix(Tset, e % Nmesh(1), spA(e % Nmesh(1)), e % Nout(1), xi)
            call addNewInterpolationMatrix(Tset, e % Nmesh(2), spA(e % Nmesh(2)), e % Nout(2), eta)   ! TODO: check why it was Nmesh(1)
            call addNewInterpolationMatrix(Tset, e % Nmesh(3), spA(e % Nmesh(3)), e % Nout(3), zeta)  ! TODO: check why it was Nmesh(1)

!
!           Construct interpolation matrices for the solution
!           -------------------------------------------------
            call addNewInterpolationMatrix(Tset, e % Nsol(1), spA(e % Nsol(1)), e % Nout(1), xi)
            call addNewInterpolationMatrix(Tset, e % Nsol(2), spA(e % Nsol(2)), e % Nout(2), eta)     ! TODO: check why it was Nsol(1)
            call addNewInterpolationMatrix(Tset, e % Nsol(3), spA(e % Nsol(3)), e % Nout(3), zeta)    ! TODO: check why it was Nsol(1)
!
!           Perform interpolation
!           ---------------------
            call ProjectStorageHomogeneousPoints(e, Tset(e % Nout(1), e % Nmesh(1)) % T, &
                                                    Tset(e % Nout(2), e % Nmesh(2)) % T, &
                                                    Tset(e % Nout(3), e % Nmesh(3)) % T, &
                                                     Tset(e % Nout(1), e % Nsol(1)) % T, &
                                                     Tset(e % Nout(2), e % Nsol(2)) % T, &
                                                     Tset(e % Nout(3), e % Nsol(3)) % T    )


            end associate
         end do
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
         write(fid,'(A,A)') 'VARIABLES = "x","y","z","Umean","Vmean","Wmean","Sxx","Syy","Szz","Sxy","Sxz","Syz"'
!
!        Write elements
!        --------------
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )

            call WriteElementToTecplot(fid, e, mesh % refs)
            end associate
         end do
!
!        Write boundaries
!        ----------------
         if (hasBoundaries) then
            do bID=1, size (mesh % boundaries)
               call WriteBoundaryToTecplot(fid, mesh % boundaries(bID), mesh % elements)
            end do
         end if
!
!        Close the file
!        --------------
         close(fid)

      end subroutine Stats2Plt_Homogeneous

      subroutine ProjectStorageHomogeneousPoints(e, TxMesh, TyMesh, TzMesh, TxSol, TySol, TzSol)
         use Storage
         use NodalStorageClass
         implicit none
         type(Element_t)     :: e
         real(kind=RP),       intent(in)  :: TxMesh(0:e % Nout(1), 0:e % Nmesh(1))
         real(kind=RP),       intent(in)  :: TyMesh(0:e % Nout(2), 0:e % Nmesh(2))
         real(kind=RP),       intent(in)  :: TzMesh(0:e % Nout(3), 0:e % Nmesh(3))
         real(kind=RP),       intent(in)  :: TxSol(0:e % Nout(1), 0:e % Nsol(1))
         real(kind=RP),       intent(in)  :: TySol(0:e % Nout(2), 0:e % Nsol(2))
         real(kind=RP),       intent(in)  :: TzSol(0:e % Nout(3), 0:e % Nsol(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j, k, iVar, l, m, n
!
!        Project mesh
!        ------------         
         allocate( e % xOut(1:3,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
         e % xOut = 0.0_RP

         do n = 0, e % Nmesh(3) ; do m = 0, e % Nmesh(2) ; do l = 0, e % Nmesh(1)
            do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
               e % xOut(:,i,j,k) = e % xOut(:,i,j,k) + e % x(:,l,m,n) * TxMesh(i,l) * TyMesh(j,m) * TzMesh(k,n)
            end do            ; end do            ; end do
         end do            ; end do            ; end do

!
!        Project the solution
!        --------------------
         allocate( e % statsout(1:9,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
         e % statsout = 0.0_RP

         do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
            do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
               e % statsout(:,i,j,k) = e % statsout(:,i,j,k) + e % stats(:,l,m,n) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
            end do            ; end do            ; end do
         end do            ; end do            ; end do

      end subroutine ProjectStorageHomogeneousPoints
!
!/////////////////////////////////////////////////////////////////////////////
!
!     Write solution
!     --------------
!
!/////////////////////////////////////////////////////////////////////////////
!
      subroutine WriteElementToTecplot(fid,e,refs)
         use Storage
         use NodalStorageClass
         use prolongMeshAndSolution
         use SolutionFile
         use StatisticsMonitor
         implicit none
         integer,            intent(in)    :: fid
         type(Element_t),    intent(inout) :: e 
         real(kind=RP),      intent(in)    :: refs(NO_OF_SAVED_REFS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: i,j,k,var
         character(len=LINE_LENGTH) :: formatout
!
!        Get output variables
!        --------------------
         allocate (e % outputVars(1:9,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) ) 
         do k = 0, e % Nout(3)   ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
            e % outputVars(1,i,j,k) = e % statsout(U,i,j,k)
            e % outputVars(2,i,j,k) = e % statsout(V,i,j,k)
            e % outputVars(3,i,j,k) = e % statsout(W,i,j,k)
            e % outputVars(4,i,j,k) = e % statsout(UU,i,j,k) - POW2(e % statsout(U,i,j,k))
            e % outputVars(5,i,j,k) = e % statsout(VV,i,j,k) - POW2(e % statsout(V,i,j,k))
            e % outputVars(6,i,j,k) = e % statsout(WW,i,j,k) - POW2(e % statsout(W,i,j,k))
            e % outputVars(7,i,j,k) = e % statsout(UV,i,j,k) - e % statsout(U,i,j,k) * e % statsout(V,i,j,k)
            e % outputVars(8,i,j,k) = e % statsout(UW,i,j,k) - e % statsout(U,i,j,k) * e % statsout(W,i,j,k)
            e % outputVars(9,i,j,k) = e % statsout(VW,i,j,k) - e % statsout(V,i,j,k) * e % statsout(W,i,j,k)
         end do                  ; end do                ; end do
!
!        Write variables
!        ---------------        
         write(fid,'(A,I0,A,I0,A,I0,A)') "ZONE I=",e % Nout(1)+1,", J=",e % Nout(2)+1, &
                                            ", K=",e % Nout(3)+1,", F=POINT"

         formatout = getFormat()

         do k = 0, e % Nout(3)   ; do j = 0, e % Nout(2)    ; do i = 0, e % Nout(1)
            write(fid,trim(formatout)) e % xOut(:,i,j,k), e % outputVars(:,i,j,k)
         end do               ; end do                ; end do

      end subroutine WriteElementToTecplot

      character(len=LINE_LENGTH) function getFormat()
         use OutputVariables
         implicit none

         getFormat = ""

         write(getFormat,'(A,I0,A,A)') "(",12,PRECISION_FORMAT,")"

      end function getFormat
      
end module Stats2PltModule