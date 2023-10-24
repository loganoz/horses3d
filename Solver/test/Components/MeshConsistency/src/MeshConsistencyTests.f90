module MeshTests
   use FTAssertions
   use SMConstants
   use FaceClass
   use ElementClass
   use MeshTypes
   use HexMeshClass
   use MeshConsistencySetup
   use ReadMeshFile
   use NodalStorageClass, only: GAUSS, NodalStorage
   use PhysicsStorage

   private

   public mesh

   public CheckVolume, CheckZoneSurfaces, CheckCoordinatesConsistency
   public CheckNormalConsistency, CheckScalConsistency, OpenNextMesh

   type(HexMesh)             :: mesh
   integer                    :: currentMesh = 0

   contains
      subroutine OpenNextMesh()
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer              :: no_of_elements
         integer, allocatable :: Nx(:), Ny(:), Nz(:)
         logical              :: success

         if ( currentMesh .ne. 0 ) then         
            call mesh % Destruct()
         end if 

         currentMesh = currentMesh + 1
         no_of_elements = NumOfElemsFromMeshFile( meshfileNames(currentMesh) )
         allocate(Nx(no_of_elements), Ny(no_of_elements), Nz(no_of_elements))
         Nx = 5
         Ny = 5
         Nz = 5

         call constructMeshFromFile( mesh, meshfileNames(currentMesh), GAUSS, Nx, Ny, Nz, .true. , 0, .false., success )

      end subroutine OpenNextMesh

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
         do eID = 1, mesh % no_of_elements
            associate(e => mesh % elements(eID))
            associate(spAxi   => NodalStorage(e % Nxyz(1)), &
                      spAeta  => NodalStorage(e % Nxyz(2)), &
                      spAzeta => NodalStorage(e % Nxyz(3)) )
            do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               volume = volume +  spAxi % w(i) * spAeta % w(j) * spAzeta % w(k) &
                                 * e % geom % jacobian(i,j,k) 
            end do                  ; end do                ; end do
            end associate
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

         do zID = 1, size(mesh % zones)
            surface = 0.0_RP
            do nF = 1, mesh % zones(zID) % no_of_faces
               fID = mesh % zones(zID) % faces(nF)
               associate(f => mesh % faces(fID))
               associate(spAxi   => NodalStorage(f % Nf(1)), &
                         spAeta  => NodalStorage(f % Nf(2)))
                  do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
                     surface = surface + spAxi % w(i) * spAeta % w(j) * f % geom % jacobian(i,j)
                  end do               ; end do
               end associate
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

         do fID = 1, size(mesh % faces)
            if ( mesh % faces(fID) % faceType .eq. HMESH_INTERIOR ) then
               eIDL = mesh % faces(fID) % elementIDs(1)
               eIDR = mesh % faces(fID) % elementIDs(2)
               call CheckElementAndFaceCoordinatesConsistency(f = mesh % faces(fID), &
                                                              eL = mesh % elements(eIDL), & 
                                                              eR = mesh % elements(eIDR),&
                                                              localError = localError)
            else
               eIDL = mesh % faces(fID) % elementIDs(1)
               call CheckElementAndFaceCoordinatesConsistency(f = mesh % faces(fID), &
                                                              eL = mesh % elements(eIDL), & 
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
         associate(spAxi   => NodalStorage(eL % Nxyz(1)), &
                   spAeta  => NodalStorage(eL % Nxyz(2)), &
                   spAzeta => NodalStorage(eL % Nxyz(3)) )

         select case( f % elementSide(1) )
         case (EFRONT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               XelL(:,i,k) = XelL(:,i,k) + eL % geom % X(:,i,j,k) * spAeta % v(j,FRONT)
            end do                  ; end do                   ; end do

         case (EBACK)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               XelL(:,i,k) = XelL(:,i,k) + eL % geom % X(:,i,j,k) * spAeta % v(j,BACK)
            end do                  ; end do                   ; end do
   
         case (EBOTTOM)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               XelL(:,i,j) = XelL(:,i,j) + eL % geom % X(:,i,j,k) * spAzeta % v(k,BOTTOM)
            end do                  ; end do                   ; end do

         case (ERIGHT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               XelL(:,j,k) = XelL(:,j,k) + eL % geom % X(:,i,j,k) * spAxi % v(i,RIGHT)
            end do                  ; end do                   ; end do

         case (ETOP)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               XelL(:,i,j) = XelL(:,i,j) + eL % geom % X(:,i,j,k) * spAzeta % v(k,TOP)
            end do                  ; end do                   ; end do

         case (ELEFT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               XelL(:,j,k) = XelL(:,j,k) + eL % geom % X(:,i,j,k) * spAxi % v(i,LEFT)
            end do                  ; end do                   ; end do

         end select
         end associate
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
         associate(spAxi   => NodalStorage(eR % Nxyz(1)), &
                   spAeta  => NodalStorage(eR % Nxyz(2)), &
                   spAzeta => NodalStorage(eR % Nxyz(3)) )
         select case( f % elementSide(2) )
         case (EFRONT)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               XelR(:,i,k) = XelR(:,i,k) + eR % geom % X(:,i,j,k) * spAeta % v(j,FRONT)
            end do                  ; end do                   ; end do

         case (EBACK)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               XelR(:,i,k) = XelR(:,i,k) + eR % geom % X(:,i,j,k) * spAeta % v(j,BACK)
            end do                  ; end do                   ; end do
   
         case (EBOTTOM)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               XelR(:,i,j) = XelR(:,i,j) + eR % geom % X(:,i,j,k) * spAzeta % v(k,BOTTOM)
            end do                  ; end do                   ; end do

         case (ERIGHT)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               XelR(:,j,k) = XelR(:,j,k) + eR % geom % X(:,i,j,k) * spAxi % v(i,RIGHT)
            end do                  ; end do                   ; end do

         case (ETOP)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               XelR(:,i,j) = XelR(:,i,j) + eR % geom % X(:,i,j,k) * spAzeta % v(k,TOP)
            end do                  ; end do                   ; end do

         case (ELEFT)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               XelR(:,j,k) = XelR(:,j,k) + eR % geom % X(:,i,j,k) * spAxi % v(i,LEFT)
            end do                  ; end do                   ; end do

         end select
         end associate

         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            call leftIndexes2Right(i,j,f % Nf(1), f % Nf(2), f % rotation, ii, jj)
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

         do fID = 1, size(mesh % faces)
            if ( mesh % faces(fID) % faceType .eq. HMESH_INTERIOR ) then
               eIDL = mesh % faces(fID) % elementIDs(1)
               eIDR = mesh % faces(fID) % elementIDs(2)
               call CheckElementAndFaceNormalConsistency(f = mesh % faces(fID), &
                                                              eL = mesh % elements(eIDL), & 
                                                              eR = mesh % elements(eIDR),&
                                                              localError = localError)
            else
               eIDL = mesh % faces(fID) % elementIDs(1)
               call CheckElementAndFaceNormalConsistency(f = mesh % faces(fID), &
                                                              eL = mesh % elements(eIDL), & 
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
         associate(spAxi   => NodalStorage(eL % Nxyz(1)), &
                   spAeta  => NodalStorage(eL % Nxyz(2)), &
                   spAzeta => NodalStorage(eL % Nxyz(3)) )
         select case( f % elementSide(1) )
         case (EFRONT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,i,k) = SelL(:,i,k) + eL % geom % jGradEta(:,i,j,k) * spAeta % v(j,FRONT)
            end do                  ; end do                   ; end do

            SelL = -SelL

         case (EBACK)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,i,k) = SelL(:,i,k) + eL % geom % jGradEta(:,i,j,k) * spAeta % v(j,BACK)
            end do                  ; end do                   ; end do
   
         case (EBOTTOM)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,i,j) = SelL(:,i,j) + eL % geom % jGradZeta(:,i,j,k) * spAzeta % v(k,BOTTOM)
            end do                  ; end do                   ; end do

            SelL = -SelL

         case (ERIGHT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,j,k) = SelL(:,j,k) + eL % geom % jGradXi(:,i,j,k) * spAxi % v(i,RIGHT)
            end do                  ; end do                   ; end do

         case (ETOP)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,i,j) = SelL(:,i,j) + eL % geom % jGradZeta(:,i,j,k) * spAzeta % v(k,TOP)
            end do                  ; end do                   ; end do

         case (ELEFT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,j,k) = SelL(:,j,k) + eL % geom % jGradXi(:,i,j,k) * spAxi % v(i,LEFT)
            end do                  ; end do                   ; end do

            SelL = -SelL

         end select
         end associate
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
         associate(spAxi   => NodalStorage(eR % Nxyz(1)), &
                   spAeta  => NodalStorage(eR % Nxyz(2)), &
                   spAzeta => NodalStorage(eR % Nxyz(3)) )
         select case( f % elementSide(2) )
         case (EFRONT)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,i,k) = SelR(:,i,k) + eR % geom % jGradEta(:,i,j,k) * spAeta % v(j,FRONT)
            end do                  ; end do                   ; end do

            SelR = -SelR

         case (EBACK)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,i,k) = SelR(:,i,k) + eR % geom % jGradEta(:,i,j,k) * spAeta % v(j,BACK)
            end do                  ; end do                   ; end do
   
         case (EBOTTOM)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,i,j) = SelR(:,i,j) + eR % geom % jGradZeta(:,i,j,k) * spAzeta % v(k,BOTTOM)
            end do                  ; end do                   ; end do

            SelR = -SelR

         case (ERIGHT)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,j,k) = SelR(:,j,k) + eR % geom % jGradXi(:,i,j,k) * spAxi % v(i,RIGHT)
            end do                  ; end do                   ; end do

         case (ETOP)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,i,j) = SelR(:,i,j) + eR % geom % jGradZeta(:,i,j,k) * spAzeta % v(k,TOP)
            end do                  ; end do                   ; end do

         case (ELEFT)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,j,k) = SelR(:,j,k) + eR % geom % jGradXi(:,i,j,k) * spAxi % v(i,LEFT)
            end do                  ; end do                   ; end do

            SelR = -SelR

         end select
         end associate
!
!        --------------------
!        Perform the rotation
!        --------------------
!
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            call leftIndexes2Right(i,j,f % Nf(1), f % Nf(2), f % rotation, ii, jj)
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

         do fID = 1, size(mesh % faces)
            if ( mesh % faces(fID) % faceType .eq. HMESH_INTERIOR ) then
               eIDL = mesh % faces(fID) % elementIDs(1)
               eIDR = mesh % faces(fID) % elementIDs(2)
               call CheckElementAndFaceScalConsistency(f = mesh % faces(fID), &
                                                              eL = mesh % elements(eIDL), & 
                                                              eR = mesh % elements(eIDR),&
                                                              localError = localError)
            else
               eIDL = mesh % faces(fID) % elementIDs(1)
               call CheckElementAndFaceScalConsistency(f = mesh % faces(fID), &
                                                              eL = mesh % elements(eIDL), & 
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
         real(kind=RP)  :: jacobianL(0:f % NelLeft(1), 0:f % NelLeft(2))
         real(kind=RP)  :: jacobianR(0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP)  :: SelR(NDIM, 0:f % NelRight(1), 0:f % NelRight(2))
         real(kind=RP)  :: SfR(NDIM, 0:f % NfRight(1), 0:f % NfRight(2))
         real(kind=RP)  :: Sf(NDIM, 0:f % Nf(1), 0:f % Nf(2))
         character(len=72)    :: msg

         SelL = 0.0_RP
         
         associate(spAxi   => NodalStorage(eL % Nxyz(1)), &
                   spAeta  => NodalStorage(eL % Nxyz(2)), &
                   spAzeta => NodalStorage(eL % Nxyz(3)) )
         select case( f % elementSide(1) )
         case (EFRONT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,i,k) = SelL(:,i,k) + eL % geom % jGradEta(:,i,j,k) * spAeta % v(j,FRONT)
            end do                  ; end do                   ; end do

            SelL = -SelL

         case (EBACK)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,i,k) = SelL(:,i,k) + eL % geom % jGradEta(:,i,j,k) * spAeta % v(j,BACK)
            end do                  ; end do                   ; end do
   
         case (EBOTTOM)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,i,j) = SelL(:,i,j) + eL % geom % jGradZeta(:,i,j,k) * spAzeta % v(k,BOTTOM)
            end do                  ; end do                   ; end do

            SelL = -SelL

         case (ERIGHT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,j,k) = SelL(:,j,k) + eL % geom % jGradXi(:,i,j,k) * spAxi % v(i,RIGHT)
            end do                  ; end do                   ; end do

         case (ETOP)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,i,j) = SelL(:,i,j) + eL % geom % jGradZeta(:,i,j,k) * spAzeta % v(k,TOP)
            end do                  ; end do                   ; end do

         case (ELEFT)
            do k = 0, eL % Nxyz(3)  ; do j = 0, eL % Nxyz(2)   ; do i = 0, eL % Nxyz(1)
               SelL(:,j,k) = SelL(:,j,k) + eL % geom % jGradXi(:,i,j,k) * spAxi % v(i,LEFT)
            end do                  ; end do                   ; end do

            SelL = -SelL

         end select
         end associate


         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            jacobianL(i,j) = norm2(SelL(:,i,j))
         end do               ; end do
!
!        -------------
!        Perform tests
!        -------------
!
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            write(msg,'(A,I0,A,I0,A,I0,A,I0,A)') "Surface Jacobian (face ",f % ID, &
                     ", element ", eL % eID," localcoords:",i,",",j,")"
            call FTAssertEqual( expectedValue = f % geom % jacobian(i,j), &
                                actualValue = jacobianL(i,j), &
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
         associate(spAxi   => NodalStorage(eR % Nxyz(1)), &
                   spAeta  => NodalStorage(eR % Nxyz(2)), &
                   spAzeta => NodalStorage(eR % Nxyz(3)) )
         select case( f % elementSide(2) )
         case (EFRONT)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,i,k) = SelR(:,i,k) + eR % geom % jGradEta(:,i,j,k) * spAeta % v(j,FRONT)
            end do                  ; end do                   ; end do

            SelR = -SelR

         case (EBACK)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,i,k) = SelR(:,i,k) + eR % geom % jGradEta(:,i,j,k) * spAeta % v(j,BACK)
            end do                  ; end do                   ; end do
   
         case (EBOTTOM)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,i,j) = SelR(:,i,j) + eR % geom % jGradZeta(:,i,j,k) * spAzeta % v(k,BOTTOM)
            end do                  ; end do                   ; end do

            SelR = -SelR

         case (ERIGHT)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,j,k) = SelR(:,j,k) + eR % geom % jGradXi(:,i,j,k) * spAxi % v(i,RIGHT)
            end do                  ; end do                   ; end do

         case (ETOP)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,i,j) = SelR(:,i,j) + eR % geom % jGradZeta(:,i,j,k) * spAzeta % v(k,TOP)
            end do                  ; end do                   ; end do

         case (ELEFT)
            do k = 0, eR % Nxyz(3)  ; do j = 0, eR % Nxyz(2)   ; do i = 0, eR % Nxyz(1)
               SelR(:,j,k) = SelR(:,j,k) + eR % geom % jGradXi(:,i,j,k) * spAxi % v(i,LEFT)
            end do                  ; end do                   ; end do

            SelR = -SelR

         end select
         end associate
!
!        --------------------
!        Perform the rotation
!        --------------------
!
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            call leftIndexes2Right(i,j,f % Nf(1), f % Nf(2), f % rotation, ii, jj)
            SfR(:,i,j) = SelR(:,ii,jj)
         end do               ; end do

         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            jacobianR(i,j) = norm2(SfR(:,i,j))
         end do               ; end do
!
!        -------------
!        Perform tests
!        -------------
!
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            write(msg,'(A,I0,A,I0,A,I0,A,I0,A)') "Surface Jacobian (face ",f % ID, &
                     ", element ", eL % eID," localcoords:",i,",",j,")"
            call FTAssertEqual( expectedValue = f % geom % jacobian(i,j), &
                                actualValue = jacobianR(i,j), &
                                tol = 1.0e-13_RP, &
                                msg = msg)
         end do               ; end do

      end subroutine CheckElementAndFaceScalConsistency
end module MeshTests