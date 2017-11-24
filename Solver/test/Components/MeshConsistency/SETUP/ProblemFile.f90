module MeshTests
   use FTAssertions
   use SMConstants
   use HexMeshClass
   use PhysicsStorage

   private

   public mesh_

   public CheckVolume, CheckZoneSurfaces, CheckCoordinatesConsistency
   public CheckNormalConsistency, CheckScalConsistency

   class(HexMesh), pointer    :: mesh_ => NULL()

   contains
      subroutine CheckVolume()
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: eID, i, j, k
         real(kind=RP)  :: volume

         volume = 0.0_RP
         do eID = 1, mesh_ % no_of_elements
            associate(e => mesh_ % elements(eID))
            do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               volume = volume +  e % spAxi % w(i) * e % spAeta % w(j) * e % spAzeta % w(k) &
                                 * e % geom % jacobian(i,j,k) 
            end do                  ; end do                ; end do
            end associate
         end do

         call FTAssertEqual(expectedValue = 125.0_RP, &
                            actualValue = volume, &
                            tol = 1.0e-13_RP, &
                            msg = "Mesh volume")

      end subroutine CheckVolume

      subroutine CheckZoneSurfaces()
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: zID, nF, fID, i, j
         real(kind=RP)     :: surface

         do zID = 1, size(mesh_ % zones)
            surface = 0.0_RP
            do nF = 1, mesh_ % zones(zID) % no_of_faces
               fID = mesh_ % zones(zID) % faces(nF)
               associate(f => mesh_ % faces(fID))
                  do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                     surface = surface + f % spAxi % w(i) * f % spAeta % w(j) * f % geom % scal(i,j)
                  end do               ; end do
               end associate
            end do

            call FTAssertEqual(expectedValue = 25.0_RP, &
                            actualValue = surface, &
                            tol = 1.0e-13_RP, &
                            msg = "Mesh volume")

         end do

      end subroutine CheckZoneSurfaces

      subroutine CheckCoordinatesConsistency()
         use FaceClass
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: eIDL, eIDR, fID
         real(kind=RP)  :: localError

         do fID = 1, size(mesh_ % faces)
            if ( mesh_ % faces(fID) % faceType .eq. HMESH_INTERIOR ) then
               eIDL = mesh_ % faces(fID) % elementIDs(1)
               eIDR = mesh_ % faces(fID) % elementIDs(2)
               call CheckElementAndFaceCoordinatesConsistency(f = mesh_ % faces(fID), &
                                                              eL = mesh_ % elements(eIDL), & 
                                                              eR = mesh_ % elements(eIDR),&
                                                              localError = localError)
            else
               eIDL = mesh_ % faces(fID) % elementIDs(1)
               call CheckElementAndFaceCoordinatesConsistency(f = mesh_ % faces(fID), &
                                                              eL = mesh_ % elements(eIDL), & 
                                                              localError = localError)

            end if
         end do


      end subroutine CheckCoordinatesConsistency

      subroutine CheckElementAndFaceCoordinatesConsistency(f, eL, eR, localError)
         use FaceClass
         use ElementClass
         use ElementConnectivityDefinitions
         implicit none
         type(Face), intent(in)     :: f
         type(Element), intent(in) :: eL
         type(Element), intent(in), optional  :: eR
         real(kind=RP), intent(out) :: localError
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, k, ii, jj
         real(kind=RP)  :: XelL(NDIM, 0:f % NelLeft(1), 0:f % NelLeft(2))
         real(kind=RP)  :: XelR(NDIM, 0:f % NelRight(1), 0:f % NelRight(2))
         real(kind=RP)  :: XfR(NDIM, 0:f % NfRight(1), 0:f % NfRight(2))
         character(len=72) :: msg
!  TODO account for p-Adaption!         
!
!        ******************
!        Check Left element
!        ******************
!
         XelL = 0.0_RP

         select case( f % elementSide(1) )
         case (EFRONT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               XelL(:,i,k) = XelL(:,i,k) + eL % geom % X(:,i,j,k) * eL % spAeta % v(j,FRONT)
            end do                  ; end do                   ; end do

         case (EBACK)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               XelL(:,i,k) = XelL(:,i,k) + eL % geom % X(:,i,j,k) * eL % spAeta % v(j,BACK)
            end do                  ; end do                   ; end do
   
         case (EBOTTOM)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               XelL(:,i,j) = XelL(:,i,j) + eL % geom % X(:,i,j,k) * eL % spAzeta % v(k,BOTTOM)
            end do                  ; end do                   ; end do

         case (ERIGHT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               XelL(:,j,k) = XelL(:,j,k) + eL % geom % X(:,i,j,k) * eL % spAxi % v(i,RIGHT)
            end do                  ; end do                   ; end do

         case (ETOP)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               XelL(:,i,j) = XelL(:,i,j) + eL % geom % X(:,i,j,k) * eL % spAzeta % v(k,TOP)
            end do                  ; end do                   ; end do

         case (ELEFT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               XelL(:,j,k) = XelL(:,j,k) + eL % geom % X(:,i,j,k) * eL % spAxi % v(i,LEFT)
            end do                  ; end do                   ; end do

         end select
!
!        -------------
!        Perform tests
!        -------------
!
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            write(msg,'(A,I0,A,I0,A,I0,A,I0,A)') "Face coordinates interpolation (face ",f % ID, &
                     ", element ", eL % eID," localcoords:",i,",",j,")"
            call FTAssertEqual( expectedValue = f % geom % x(1,i,j), &
                                actualValue = XelL(1,i,j), &
                                tol = 1.0e-13_RP, &
                                msg = msg)
            call FTAssertEqual( expectedValue = f % geom % x(2,i,j), &
                                actualValue = XelL(2,i,j), &
                                tol = 1.0e-13_RP, &
                                msg = msg)
            call FTAssertEqual( expectedValue = f % geom % x(3,i,j), &
                                actualValue = XelL(3,i,j), &
                                tol = 1.0e-13_RP, &
                                msg = msg)
         end do               ; end do
   


!
!        *******************
!        Check Right element
!        *******************
!
         if ( .not. present(eR) ) return

         XelR = 0.0_RP

         select case( f % elementSide(2) )
         case (EFRONT)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               XelR(:,i,k) = XelR(:,i,k) + eR % geom % X(:,i,j,k) * eR % spAeta % v(j,FRONT)
            end do                  ; end do                   ; end do

         case (EBACK)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               XelR(:,i,k) = XelR(:,i,k) + eR % geom % X(:,i,j,k) * eR % spAeta % v(j,BACK)
            end do                  ; end do                   ; end do
   
         case (EBOTTOM)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               XelR(:,i,j) = XelR(:,i,j) + eR % geom % X(:,i,j,k) * eR % spAzeta % v(k,BOTTOM)
            end do                  ; end do                   ; end do

         case (ERIGHT)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               XelR(:,j,k) = XelR(:,j,k) + eR % geom % X(:,i,j,k) * eR % spAxi % v(i,RIGHT)
            end do                  ; end do                   ; end do

         case (ETOP)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               XelR(:,i,j) = XelR(:,i,j) + eR % geom % X(:,i,j,k) * eR % spAzeta % v(k,TOP)
            end do                  ; end do                   ; end do

         case (ELEFT)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               XelR(:,j,k) = XelR(:,j,k) + eR % geom % X(:,i,j,k) * eR % spAxi % v(i,LEFT)
            end do                  ; end do                   ; end do

         end select

         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            call iijjIndexes(i,j,f % Nf(1), f % Nf(2), f % rotation, ii, jj)
            XfR(:,i,j) = XelR(:,ii,jj)
         end do               ; end do
!
!        -------------
!        Perform tests
!        -------------
!
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            write(msg,'(A,I0,A,I0,A,I0,A,I0,A)') "Face coordinates interpolation (face ",f % ID, &
                     ", element ", eR % eID," localcoords:",i,",",j,")"
            call FTAssertEqual( expectedValue = f % geom % x(1,i,j), &
                                actualValue = XfR(1,i,j), &
                                tol = 1.0e-13_RP, &
                                msg = msg)
            call FTAssertEqual( expectedValue = f % geom % x(2,i,j), &
                                actualValue = XfR(2,i,j), &
                                tol = 1.0e-13_RP, &
                                msg = msg)
            call FTAssertEqual( expectedValue = f % geom % x(3,i,j), &
                                actualValue = XfR(3,i,j), &
                                tol = 1.0e-13_RP, &
                                msg = msg)
         end do               ; end do


      end subroutine CheckElementAndFaceCoordinatesConsistency

      subroutine CheckNormalConsistency()
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: eIDL, eIDR, fID
         real(kind=RP)  :: localError

         do fID = 1, size(mesh_ % faces)
            if ( mesh_ % faces(fID) % faceType .eq. HMESH_INTERIOR ) then
               eIDL = mesh_ % faces(fID) % elementIDs(1)
               eIDR = mesh_ % faces(fID) % elementIDs(2)
               call CheckElementAndFaceNormalConsistency(f = mesh_ % faces(fID), &
                                                              eL = mesh_ % elements(eIDL), & 
                                                              eR = mesh_ % elements(eIDR),&
                                                              localError = localError)
            else
               eIDL = mesh_ % faces(fID) % elementIDs(1)
               call CheckElementAndFaceNormalConsistency(f = mesh_ % faces(fID), &
                                                              eL = mesh_ % elements(eIDL), & 
                                                              localError = localError)

            end if
         end do



      end subroutine CheckNormalConsistency

      subroutine CheckElementAndFaceNormalConsistency(f, eL, eR, localError)
         implicit none
         type(Face), intent(in)              :: f
         type(Element), intent(in)           :: eL
         type(Element), intent(in), optional :: eR
         real(kind=RP), intent(out)          :: localError
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, k, ii, jj
         real(kind=RP)  :: SelL(NDIM, 0:f % NelLeft(1), 0:f % NelLeft(2))
         real(kind=RP)  :: SelR(NDIM, 0:f % NelRight(1), 0:f % NelRight(2))
         real(kind=RP)  :: SfR(NDIM, 0:f % NfRight(1), 0:f % NfRight(2))
         character(len=72) :: msg

         SelL = 0.0_RP

         select case( f % elementSide(1) )
         case (EFRONT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,i,k) = SelL(:,i,k) + eL % geom % jGradEta(:,i,j,k) * eL % spAeta % v(j,FRONT)
            end do                  ; end do                   ; end do

            SelL = -SelL

         case (EBACK)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,i,k) = SelL(:,i,k) + eL % geom % jGradEta(:,i,j,k) * eL % spAeta % v(j,BACK)
            end do                  ; end do                   ; end do
   
         case (EBOTTOM)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,i,j) = SelL(:,i,j) + eL % geom % jGradZeta(:,i,j,k) * eL % spAzeta % v(k,BOTTOM)
            end do                  ; end do                   ; end do

            SelL = -SelL

         case (ERIGHT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,j,k) = SelL(:,j,k) + eL % geom % jGradXi(:,i,j,k) * eL % spAxi % v(i,RIGHT)
            end do                  ; end do                   ; end do

         case (ETOP)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,i,j) = SelL(:,i,j) + eL % geom % jGradZeta(:,i,j,k) * eL % spAzeta % v(k,TOP)
            end do                  ; end do                   ; end do

         case (ELEFT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,j,k) = SelL(:,j,k) + eL % geom % jGradXi(:,i,j,k) * eL % spAxi % v(i,LEFT)
            end do                  ; end do                   ; end do

            SelL = -SelL

         end select
!
!        ---------------------
!        Get the normal vector
!        ---------------------
!
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            SelL(:,i,j) = SelL(:,i,j) / norm2(SelL(:,i,j))
         end do               ; end do
!
!        -------------
!        Perform tests
!        -------------
!
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            write(msg,'(A,I0,A,I0,A,I0,A,I0,A)') "Face normals (face ",f % ID, &
                     ", element ", eL % eID," localcoords:",i,",",j,")"
            call FTAssertEqual( expectedValue = f % geom % normal(1,i,j), &
                                actualValue = SelL(1,i,j), &
                                tol = 1.0e-13_RP, &
                                msg = msg)
            call FTAssertEqual( expectedValue = f % geom % normal(2,i,j), &
                                actualValue = SelL(2,i,j), &
                                tol = 1.0e-13_RP, &
                                msg = msg)
            call FTAssertEqual( expectedValue = f % geom % normal(3,i,j), &
                                actualValue = SelL(3,i,j), &
                                tol = 1.0e-13_RP, &
                                msg = msg)
         end do               ; end do
!
!        ***********************
!        Check the right element
!        ***********************
!
         if ( .not. present(eR) ) return

         SelR = 0.0_RP

         select case( f % elementSide(2) )
         case (EFRONT)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,i,k) = SelR(:,i,k) + eR % geom % jGradEta(:,i,j,k) * eR % spAeta % v(j,FRONT)
            end do                  ; end do                   ; end do

            SelR = -SelR

         case (EBACK)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,i,k) = SelR(:,i,k) + eR % geom % jGradEta(:,i,j,k) * eR % spAeta % v(j,BACK)
            end do                  ; end do                   ; end do
   
         case (EBOTTOM)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,i,j) = SelR(:,i,j) + eR % geom % jGradZeta(:,i,j,k) * eR % spAzeta % v(k,BOTTOM)
            end do                  ; end do                   ; end do

            SelR = -SelR

         case (ERIGHT)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,j,k) = SelR(:,j,k) + eR % geom % jGradXi(:,i,j,k) * eR % spAxi % v(i,RIGHT)
            end do                  ; end do                   ; end do

         case (ETOP)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,i,j) = SelR(:,i,j) + eR % geom % jGradZeta(:,i,j,k) * eR % spAzeta % v(k,TOP)
            end do                  ; end do                   ; end do

         case (ELEFT)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,j,k) = SelR(:,j,k) + eR % geom % jGradXi(:,i,j,k) * eR % spAxi % v(i,LEFT)
            end do                  ; end do                   ; end do

            SelR = -SelR

         end select
!
!        --------------------
!        Perform the rotation
!        --------------------
!
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            call iijjIndexes(i,j,f % Nf(1), f % Nf(2), f % rotation, ii, jj)
            SfR(:,i,j) = SelR(:,ii,jj)
         end do               ; end do

         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            SfR(:,i,j) = SfR(:,i,j) / norm2(SfR(:,i,j))
         end do               ; end do
!
!        -------------
!        Perform tests
!        -------------
!
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            write(msg,'(A,I0,A,I0,A,I0,A,I0,A)') "Face normals (face ",f % ID, &
                     ", element ", eR % eID," localcoords:",i,",",j,")"
            call FTAssertEqual( expectedValue = f % geom % normal(1,i,j), &
                                actualValue = -SfR(1,i,j), &
                                tol = 1.0e-13_RP, &
                                msg = msg)
            call FTAssertEqual( expectedValue = f % geom % normal(2,i,j), &
                                actualValue = -SfR(2,i,j), &
                                tol = 1.0e-13_RP, &
                                msg = msg)
            call FTAssertEqual( expectedValue = f % geom % normal(3,i,j), &
                                actualValue = -SfR(3,i,j), &
                                tol = 1.0e-13_RP, &
                                msg = msg)
         end do               ; end do





      end subroutine CheckElementAndFaceNormalConsistency

      subroutine CheckScalConsistency()
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: eIDL, eIDR, fID
         real(kind=RP)  :: localError, error

         do fID = 1, size(mesh_ % faces)
            if ( mesh_ % faces(fID) % faceType .eq. HMESH_INTERIOR ) then
               eIDL = mesh_ % faces(fID) % elementIDs(1)
               eIDR = mesh_ % faces(fID) % elementIDs(2)
               call CheckElementAndFaceScalConsistency(f = mesh_ % faces(fID), &
                                                              eL = mesh_ % elements(eIDL), & 
                                                              eR = mesh_ % elements(eIDR),&
                                                              localError = localError)
            else
               eIDL = mesh_ % faces(fID) % elementIDs(1)
               call CheckElementAndFaceScalConsistency(f = mesh_ % faces(fID), &
                                                              eL = mesh_ % elements(eIDL), & 
                                                              localError = localError)

            end if
         end do

      end subroutine CheckScalConsistency

      subroutine CheckElementAndFaceScalConsistency(f, eL, eR, localError)
         implicit none
         type(Face), intent(in)              :: f
         type(Element), intent(in)           :: eL
         type(Element), intent(in), optional :: eR
         real(kind=RP), intent(out)          :: localError
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, k, ii, jj
         real(kind=RP)  :: SelL(NDIM, 0:f % NelLeft(1), 0:f % NelLeft(2))
         real(kind=RP)  :: scalL(0:f % NelLeft(1), 0:f % NelLeft(2))
         real(kind=RP)  :: SelR(NDIM, 0:f % NelRight(1), 0:f % NelRight(2))
         real(kind=RP)  :: SfR(NDIM, 0:f % NfRight(1), 0:f % NfRight(2))
         real(kind=RP)  :: Sf(NDIM, 0:f % Nf(1), 0:f % Nf(2))
         character(len=72)    :: msg

         SelL = 0.0_RP

         select case( f % elementSide(1) )
         case (EFRONT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,i,k) = SelL(:,i,k) + eL % geom % jGradEta(:,i,j,k) * eL % spAeta % v(j,FRONT)
            end do                  ; end do                   ; end do

            SelL = -SelL

         case (EBACK)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,i,k) = SelL(:,i,k) + eL % geom % jGradEta(:,i,j,k) * eL % spAeta % v(j,BACK)
            end do                  ; end do                   ; end do
   
         case (EBOTTOM)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,i,j) = SelL(:,i,j) + eL % geom % jGradZeta(:,i,j,k) * eL % spAzeta % v(k,BOTTOM)
            end do                  ; end do                   ; end do

            SelL = -SelL

         case (ERIGHT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,j,k) = SelL(:,j,k) + eL % geom % jGradXi(:,i,j,k) * eL % spAxi % v(i,RIGHT)
            end do                  ; end do                   ; end do

         case (ETOP)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,i,j) = SelL(:,i,j) + eL % geom % jGradZeta(:,i,j,k) * eL % spAzeta % v(k,TOP)
            end do                  ; end do                   ; end do

         case (ELEFT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,j,k) = SelL(:,j,k) + eL % geom % jGradXi(:,i,j,k) * eL % spAxi % v(i,LEFT)
            end do                  ; end do                   ; end do

            SelL = -SelL

         end select


         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            scalL(i,j) = norm2(SelL(:,i,j))
         end do               ; end do
!
!        -------------
!        Perform tests
!        -------------
!
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            write(msg,'(A,I0,A,I0,A,I0,A,I0,A)') "Face normals (face ",f % ID, &
                     ", element ", eL % eID," localcoords:",i,",",j,")"
            call FTAssertEqual( expectedValue = f % geom % scal(i,j), &
                                actualValue = scalL(i,j), &
                                tol = 1.0e-13_RP, &
                                msg = msg)
         end do               ; end do



      end subroutine CheckElementAndFaceScalConsistency

      
end module MeshTests
!
!////////////////////////////////////////////////////////////////////////
!
!      ProblemFile.f90
!      Created: June 26, 2015 at 8:47 AM 
!      By: David Kopriva  
!
!      The Problem File contains user defined procedures
!      that are used to "personalize" i.e. define a specific
!      problem to be solved. These procedures include initial conditions,
!      exact solutions (e.g. for tests), etc. and allow modifications 
!      without having to modify the main code.
!
!      The procedures, *even if empty* that must be defined are
!
!      UserDefinedSetUp
!      UserDefinedInitialCondition(mesh)
!      UserDefinedPeriodicOperation(mesh)
!      UserDefinedFinalize(mesh)
!      UserDefinedTermination
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedStartup
!
!        --------------------------------
!        Called before any other routines
!        --------------------------------
!
            IMPLICIT NONE  
         END SUBROUTINE UserDefinedStartup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalSetup(mesh , thermodynamics_, &
                                                 dimensionless_, &
                                                     refValues_ )
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            use TestSuiteManagerClass
            USE HexMeshClass
            use PhysicsStorage
            use MeshTests
            IMPLICIT NONE
            CLASS(HexMesh), target              :: mesh
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
            type(TestSuiteManager)              :: testSuite
            integer                             :: numberOfFailures
!
!           *************************
!           Perform tests on the mesh
!           *************************
!
            write(STD_OUT,'(/,/,A,/)') "  * Performing tests in the constructed mesh"
   
            call testSuite % init()
!
!           Set the test mesh target
!           ------------------------
            mesh_ => mesh

            call testSuite % addTestSubroutineWithName(CheckVolume,"Mesh volume")
            call testSuite % addTestSubroutineWithName(CheckZoneSurfaces,"Zone surfaces")
            call testSuite % addTestSubroutineWithName(CheckCoordinatesConsistency,&
                        "Mapping coordinates interpolation consistency")
            call testSuite % addTestSubroutineWithName(CheckNormalConsistency,&
                        "Consistency in the normal vectors")
            call testSuite % addTestSubroutineWithName(CheckScalConsistency,&
                        "Consistency in the surface Jacobians")

            call testSuite % PerformTests(numberOfFailures)
            
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedInitialCondition(mesh, thermodynamics_, &
                                                      dimensionless_, &
                                                          refValues_  )
!
!           ------------------------------------------------
!           Called to set the initial condition for the flow
!              - By default it sets an uniform initial
!                 condition.
!           ------------------------------------------------
!
            USE SMConstants
            use PhysicsStorage
            use HexMeshClass
            implicit none
            class(HexMesh)                      :: mesh
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
!
!           ---------------
!           Local variables
!           ---------------
!
            integer        :: eID, i, j, k
            real(kind=RP)  :: qq, u, v, w, p
            real(kind=RP)  :: Q(N_EQN), phi, theta

            associate ( gammaM2 => dimensionless_ % gammaM2, &
                        gamma => thermodynamics_ % gamma )
            theta = refValues_ % AOATheta*(PI/180.0_RP)
            phi   = refValues_ % AOAPhi*(PI/180.0_RP)
      
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          Ny => mesh % elements(eID) % Nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 
                  qq = 1.0_RP
                  u  = qq*cos(theta)*COS(phi)
                  v  = qq*sin(theta)*COS(phi)
                  w  = qq*SIN(phi)
      
                  Q(1) = 1.0_RP
                  p    = 1.0_RP/(gammaM2)
                  Q(2) = Q(1)*u
                  Q(3) = Q(1)*v
                  Q(4) = Q(1)*w
                  Q(5) = p/(gamma - 1._RP) + 0.5_RP*Q(1)*(u**2 + v**2 + w**2)

                  mesh % elements(eID) % storage % Q(:,i,j,k) = Q 
               end do;        end do;        end do
               end associate
            end do

            end associate
            
         END SUBROUTINE UserDefinedInitialCondition

         subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
!
!           -------------------------------------------------
!           Used to define an user defined boundary condition
!           -------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: Q(N_EQN)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
         end subroutine UserDefinedState1

         subroutine UserDefinedNeumann(x, t, nHat, U_x, U_y, U_z)
!
!           --------------------------------------------------------
!           Used to define a Neumann user defined boundary condition
!           --------------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: U_x(N_GRAD_EQN)
            real(kind=RP), intent(inout)  :: U_y(N_GRAD_EQN)
            real(kind=RP), intent(inout)  :: U_z(N_GRAD_EQN)
         end subroutine UserDefinedNeumann

!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, Monitors)
!
!           ----------------------------------------------------------
!           Called at the output interval to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            USE HexMeshClass
            use MonitorsClass
            IMPLICIT NONE
            CLASS(HexMesh)               :: mesh
            REAL(KIND=RP)                :: time
            type(Monitor_t), intent(in) :: monitors
            
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
         subroutine UserDefinedSourceTerm(mesh, time, thermodynamics_, dimensionless_, refValues_)
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
            USE HexMeshClass
            use PhysicsStorage
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            type(Thermodynamics_t),    intent(in) :: thermodynamics_
            type(Dimensionless_t),     intent(in) :: dimensionless_
            type(RefValues_t),         intent(in) :: refValues_
!
!           ---------------
!           Local variables
!           ---------------
!
            integer  :: i, j, k, eID
!
!           Usage example (by default no source terms are added)
!           ----------------------------------------------------
!           do eID = 1, mesh % no_of_elements
!              associate ( e => mesh % elements(eID) )
!              do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
!                 e % QDot(:,i,j,k) = e % QDot(:,i,j,k) + Source(:)
!              end do                  ; end do                ; end do
!           end do
   
         end subroutine UserDefinedSourceTerm
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual, thermodynamics_, &
                                                    dimensionless_, &
                                                        refValues_, &  
                                                          monitors, &
                                                       elapsedTime, &
                                                           CPUTime   )
!
!           --------------------------------------------------------
!           Called after the solution computed to allow, for example
!           error tests to be performed
!           --------------------------------------------------------
!
            USE HexMeshClass
            use PhysicsStorage
            use MonitorsClass
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
            type(Thermodynamics_t),    intent(in) :: thermodynamics_
            type(Dimensionless_t),     intent(in) :: dimensionless_
            type(RefValues_t),         intent(in) :: refValues_
            type(Monitor_t),          intent(in) :: monitors
            real(kind=RP),             intent(in)  :: elapsedTime
            real(kind=RP),             intent(in)  :: CPUTime

         END SUBROUTINE UserDefinedFinalize
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedTermination
!
!        -----------------------------------------------
!        Called at the the end of the main driver after 
!        everything else is done.
!        -----------------------------------------------
!
         IMPLICIT NONE  
      END SUBROUTINE UserDefinedTermination
      
