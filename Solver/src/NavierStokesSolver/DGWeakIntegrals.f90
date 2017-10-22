module DGWeakIntegrals
   use SMConstants
   use ElementClass
   use NodalStorageClass
   use PhysicsStorage
   use Physics
   implicit none


   private

   
   type  ScalarWeakIntegrals_t
      contains
         procedure, nopass    :: StdVolumeGreen  => ScalarWeakIntegrals_StdVolumeGreen
         procedure, nopass    :: StdFace => ScalarWeakIntegrals_StdFace
   end type ScalarWeakIntegrals_t

   type  VectorWeakIntegrals_t
      contains
         procedure, nopass    :: StdVolumeGreen  => VectorWeakIntegrals_StdVolumeGreen
         procedure, nopass    :: StdFace => VectorWeakIntegrals_StdFace
   end type VectorWeakIntegrals_t

   public ScalarWeakIntegrals_t, VectorWeakIntegrals_t
   public ScalarWeakIntegrals  , VectorWeakIntegrals
   
   type(ScalarWeakIntegrals_t)      :: ScalarWeakIntegrals
   type(VectorWeakIntegrals_t)      :: VectorWeakIntegrals
!
!  ========
   contains
!  ========
!
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!>       Weak integrals with scalar test function
!        ----------------------------------------
!/////////////////////////////////////////////////////////////////////////////////////////////
!
      function ScalarWeakIntegrals_StdVolumeGreen( e, spA, F ) result ( volInt )
         use MatrixOperations
         implicit none
         class(Element),      intent(in)  :: e
         class(NodalStorage), intent(in)  :: spA
         real(kind=RP),       intent(in)  :: F(1:NCONS, 0:spA % Nx, 0:spA % Ny, 0:spA % Nz, 1:NDIM )
         real(kind=RP)                    :: volInt(1:NCONS, 0:spA % Nx, 0:spA % Ny, 0:spA % Nz)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j, k, l

         volInt = 0.0_RP

         do k = 0, spA % Nz ; do j = 0, spA % Ny   ; do l = 0, spA % Nx ; do i = 0, spA % Nx
            volInt(:,i,j,k) = volInt(:,i,j,k) + spA % hatDx(i,l) * F(:,l,j,k,IX)
         end do             ; end do               ; end do             ; end do

         do k = 0, spA % Nz ; do l = 0, spA % Ny ; do j = 0, spA % Ny   ; do i = 0, spA % Nx
            volInt(:,i,j,k) = volInt(:,i,j,k) + spA % hatDy(j,l) * F(:,i,l,k,IY)
         end do             ; end do               ; end do             ; end do

         do l = 0, spA % Nz ; do k = 0, spA % Nz ; do j = 0, spA % Ny   ; do i = 0, spA % Nx
            volInt(:,i,j,k) = volInt(:,i,j,k) + spA % hatDz(k,l) * F(:,i,j,l,IZ)
         end do             ; end do             ; end do               ; end do

      end function ScalarWeakIntegrals_StdVolumeGreen
!
!///////////////////////////////////////////////////////////////
!
      pure function ScalarWeakIntegrals_StdFace( e, spA, F ) result ( faceInt )
         use MatrixOperations
         implicit none
         class(Element),      intent(in)     :: e
         class(NodalStorage), intent(in)     :: spA 
         real(kind=RP),       intent(in)     :: F(1:,0:,0:,1:)
         real(kind=RP)                       :: faceInt(1:NCONS, 0:spA % Nx,0:spA % Ny,0:spA % Nz)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: iXi, iEta, iZeta, iVar
!
!        ----------------
!>       Xi-contributions
!        ----------------
!
         do iZeta = 0, spA % Nz ; do iEta = 0, spA % Ny ; do iXi = 0, spa % Nx
            faceInt(:,iXi,iEta,iZeta) = F(:, iEta, iZeta,ELEFT) * spA % bx(iXi, LEFT)
         end do                 ; end do                ; end do
      
         do iZeta = 0, spA % Nz ; do iEta = 0, spA % Ny ; do iXi = 0, spa % Nx
            faceInt(:,iXi,iEta,iZeta) = faceInt(:,iXi,iEta,iZeta) + F(:, iEta, iZeta, ERIGHT) * spA % bx(iXi, RIGHT)
         end do                 ; end do                ; end do
!
!        -----------------
!>       Eta-contributions
!        -----------------
!
         do iZeta = 0, spA % Nz ; do iEta = 0, spA % Ny ; do iXi = 0, spa % Nx
            faceInt(:,iXi,iEta,iZeta) = faceInt(:,iXi,iEta,iZeta) + F(:, iXi, iZeta,EFRONT) * spA % by(iEta, LEFT)
         end do                 ; end do                ; end do

         do iZeta = 0, spA % Nz ; do iEta = 0, spA % Ny ; do iXi = 0, spa % Nx
            faceInt(:,iXi,iEta,iZeta) = faceInt(:,iXi,iEta,iZeta) + F(:, iXi, iZeta, EBACK) * spA % by(iEta, RIGHT)
         end do                 ; end do                ; end do
!
!        ------------------
!>       Zeta-contributions
!        ------------------
!
         do iZeta = 0, spA % Nz ; do iEta = 0, spA % Ny ; do iXi = 0, spa % Nx
            faceInt(:,iXi,iEta,iZeta) = faceInt(:,iXi,iEta,iZeta) + F(:, iXi, iEta, EBOTTOM) * spA % bz(iZeta, LEFT)
         end do                 ; end do                ; end do

         do iZeta = 0, spA % Nz ; do iEta = 0, spA % Ny ; do iXi = 0, spa % Nx
            faceInt(:,iXi,iEta,iZeta) = faceInt(:,iXi,iEta,iZeta) + F(:, iXi, iEta, ETOP) * spA % bz(iZeta, RIGHT)
         end do                 ; end do                ; end do

      end function ScalarWeakIntegrals_StdFace
!
!/////////////////////////////////////////////////////////////////////////////
!
!>       Weak integrals with vector test function
!        ----------------------------------------
!/////////////////////////////////////////////////////////////////////////////
!
      subroutine VectorWeakIntegrals_StdVolumeGreen( e, spA, U, volInt_x, volInt_y, volInt_z )
         use ElementClass
         use NodalStorageClass
         use Physics
         use PhysicsStorage
         implicit none
         class(Element),      intent(in)  :: e
         class(NodalStorage), intent(in)  :: spA
         real(kind=RP),       intent(in)  :: U        (N_GRAD_EQN,0:spA % Nx, 0:spA % Ny, 0:spA % Nz)
         real(kind=RP),       intent(out) :: volInt_x (N_GRAD_EQN,0:spA % Nx, 0:spA % Ny, 0:spA % Nz)
         real(kind=RP),       intent(out) :: volInt_y (N_GRAD_EQN,0:spA % Nx, 0:spA % Ny, 0:spA % Nz)
         real(kind=RP),       intent(out) :: volInt_z (N_GRAD_EQN,0:spA % Nx, 0:spA % Ny, 0:spA % Nz)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i,j,k,l
         real(kind=RP)  :: contravariantU(N_GRAD_EQN,0:spA % Nx, 0:spA % Ny, 0:spA % Nz,NDIM )

         volInt_x = 0.0_RP
         volInt_y = 0.0_RP
         volInt_z = 0.0_RP

            do k = 0, spA % Nz   ; do j = 0, spA % Ny    ; do l = 0, spA % Nx ; do i = 0, spA % Nx
               volInt_x(:,i,j,k) = volInt_x(:,i,j,k) + spA % hatDx(i,l) * e % geom % jGradXi(IX,l,j,k) * U(:,l,j,k)
               volInt_y(:,i,j,k) = volInt_y(:,i,j,k) + spA % hatDx(i,l) * e % geom % jGradXi(IY,l,j,k) * U(:,l,j,k)
               volInt_z(:,i,j,k) = volInt_z(:,i,j,k) + spA % hatDx(i,l) * e % geom % jGradXi(IZ,l,j,k) * U(:,l,j,k)
            end do               ; end do                ; end do             ; end do

            do k = 0, spA % Nz   ; do l = 0, spA % Ny ; do j = 0, spA % Ny    ; do i = 0, spA % Nx
               volInt_x(:,i,j,k) = volInt_x(:,i,j,k) + spA % hatDy(j,l) * e % geom % jGradEta(IX,i,l,k) * U(:,i,l,k)
               volInt_y(:,i,j,k) = volInt_y(:,i,j,k) + spA % hatDy(j,l) * e % geom % jGradEta(IY,i,l,k) * U(:,i,l,k)
               volInt_z(:,i,j,k) = volInt_z(:,i,j,k) + spA % hatDy(j,l) * e % geom % jGradEta(IZ,i,l,k) * U(:,i,l,k)
            end do               ; end do                ; end do    ; end do

            do l = 0, spA % Nz ; do k = 0, spA % Nz   ; do j = 0, spA % Ny    ; do i = 0, spA % Nx
               volInt_x(:,i,j,k) = volInt_x(:,i,j,k) + spA % hatDz(k,l) * e % geom % jGradZeta(IX,i,j,l) * U(:,i,j,l)
               volInt_y(:,i,j,k) = volInt_y(:,i,j,k) + spA % hatDz(k,l) * e % geom % jGradZeta(IY,i,j,l) * U(:,i,j,l)
               volInt_z(:,i,j,k) = volInt_z(:,i,j,k) + spA % hatDz(k,l) * e % geom % jGradZeta(IZ,i,j,l) * U(:,i,j,l)
            end do             ; end do               ; end do                ; end do
   
      end subroutine VectorWeakIntegrals_StdVolumeGreen
!
!/////////////////////////////////////////////////////////////////////////////////
!
      subroutine VectorWeakIntegrals_StdFace( e, spA, U, faceInt_x, faceInt_y, faceInt_z )
         use ElementClass
         use NodalStorageClass
         use Physics
         use PhysicsStorage
         implicit none
         class(Element),      intent(in)  :: e
         class(NodalStorage), intent(in)  :: spA
         real(kind=RP),       intent(in)  :: U( 1:, 0:, 0:, 1: )                                   
         real(kind=RP),       intent(out) :: faceInt_x(N_GRAD_EQN, 0:spA % Nx, 0:spA % Ny, 0:spA % Nz )  
         real(kind=RP),       intent(out) :: faceInt_y(N_GRAD_EQN, 0:spA % Nx, 0:spA % Ny, 0:spA % Nz )  
         real(kind=RP),       intent(out) :: faceInt_z(N_GRAD_EQN, 0:spA % Nx, 0:spA % Ny, 0:spA % Nz )
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: iVar, iXi, iEta, iZeta
!
!     ----------------
! >   Xi-contributions
!     ----------------
!
         do iZeta = 0, spA % Nz ; do iEta = 0, spA % Ny ; do iXi = 0, spa % Nx
            faceInt_x(:,iXi,iEta,iZeta) =   U(:, iEta, iZeta, ELEFT) * spA % bx(iXi, LEFT) &
                                             * e % geom % normal(IX,iEta,iZeta,ELEFT) * e % geom % scal(iEta,iZeta,ELEFT)

            faceInt_y(:,iXi,iEta,iZeta) =   U(:, iEta, iZeta, ELEFT) * spA % bx(iXi, LEFT) &
                                             * e % geom % normal(IY,iEta,iZeta,ELEFT) * e % geom % scal(iEta,iZeta,ELEFT)

            faceInt_z(:,iXi,iEta,iZeta) =   U(:, iEta, iZeta, ELEFT) * spA % bx(iXi, LEFT) &
                                             * e % geom % normal(IZ,iEta,iZeta,ELEFT) * e % geom % scal(iEta,iZeta,ELEFT)
         end do                 ; end do                ; end do
      

         do iZeta = 0, spA % Nz ; do iEta = 0, spA % Ny ; do iXi = 0, spa % Nx
            faceInt_x(:,iXi,iEta,iZeta) =   faceInt_x(:,iXi,iEta,iZeta) + U(:, iEta, iZeta, ERIGHT) * spA % bx(iXi, RIGHT) &
                                             * e % geom % normal(IX,iEta,iZeta,ERIGHT) * e % geom % scal(iEta,iZeta,ERIGHT)

            faceInt_y(:,iXi,iEta,iZeta) =   faceInt_y(:,iXi,iEta,iZeta) + U(:, iEta, iZeta, ERIGHT) * spA % bx(iXi, RIGHT) &
                                             * e % geom % normal(IY,iEta,iZeta,ERIGHT) * e % geom % scal(iEta,iZeta,ERIGHT)

            faceInt_z(:,iXi,iEta,iZeta) =   faceInt_z(:,iXi,iEta,iZeta) + U(:, iEta, iZeta, ERIGHT) * spA % bx(iXi, RIGHT) &
                                             * e % geom % normal(IZ,iEta,iZeta,ERIGHT) * e % geom % scal(iEta,iZeta,ERIGHT)
         end do                 ; end do                ; end do
!
!     -----------------
! >   Eta-contributions
!     -----------------
!
         do iZeta = 0, spA % Nz ; do iEta = 0, spA % Ny ; do iXi = 0, spa % Nx
            faceInt_x(:,iXi,iEta,iZeta) =   faceInt_x(:,iXi,iEta,iZeta) + U(:, iXi, iZeta, EFRONT) * spA % by(iEta, LEFT) &
                                             * e % geom % normal(IX,iXi,iZeta,EFRONT) * e % geom % scal(iXi,iZeta,EFRONT)

            faceInt_y(:,iXi,iEta,iZeta) =   faceInt_y(:,iXi,iEta,iZeta) + U(:, iXi, iZeta, EFRONT) * spA % by(iEta, LEFT) &
                                             * e % geom % normal(IY,iXi,iZeta,EFRONT) * e % geom % scal(iXi,iZeta,EFRONT)

            faceInt_z(:,iXi,iEta,iZeta) =   faceInt_z(:,iXi,iEta,iZeta) + U(:, iXi, iZeta, EFRONT) * spA % by(iEta, LEFT) &
                                             * e % geom % normal(IZ,iXi,iZeta,EFRONT) * e % geom % scal(iXi,iZeta,EFRONT)
         end do                 ; end do                ; end do

         do iZeta = 0, spA % Nz ; do iEta = 0, spA % Ny ; do iXi = 0, spa % Nx
            faceInt_x(:,iXi,iEta,iZeta) =   faceInt_x(:,iXi,iEta,iZeta) + U(:, iXi, iZeta, EBACK) * spA % by(iEta, RIGHT) &
                                             * e % geom % normal(IX,iXi,iZeta,EBACK) * e % geom % scal(iXi,iZeta,EBACK)

            faceInt_y(:,iXi,iEta,iZeta) =   faceInt_y(:,iXi,iEta,iZeta) + U(:, iXi, iZeta, EBACK) * spA % by(iEta, RIGHT) &
                                             * e % geom % normal(IY,iXi,iZeta,EBACK) * e % geom % scal(iXi,iZeta,EBACK)

            faceInt_z(:,iXi,iEta,iZeta) =   faceInt_z(:,iXi,iEta,iZeta) + U(:, iXi, iZeta, EBACK) * spA % by(iEta, RIGHT) &
                                             * e % geom % normal(IZ,iXi,iZeta,EBACK) * e % geom % scal(iXi,iZeta,EBACK)
         end do                 ; end do                ; end do
!
!     ------------------
!  >  Zeta-contributions
!     ------------------
!
         do iZeta = 0, spA % Nz ; do iEta = 0, spA % Ny ; do iXi = 0, spa % Nx
            faceInt_x(:,iXi,iEta,iZeta) =   faceInt_x(:,iXi,iEta,iZeta) + U(:, iXi, iEta, EBOTTOM) * spA % bz(iZeta, LEFT) &
                                             * e % geom % normal(IX,iXi,iEta,EBOTTOM) * e % geom % scal(iXi,iEta,EBOTTOM)

            faceInt_y(:,iXi,iEta,iZeta) =   faceInt_y(:,iXi,iEta,iZeta) + U(:, iXi, iEta, EBOTTOM) * spA % bz(iZeta, LEFT) &
                                             * e % geom % normal(IY,iXi,iEta,EBOTTOM) * e % geom % scal(iXi,iEta,EBOTTOM)

            faceInt_z(:,iXi,iEta,iZeta) =   faceInt_z(:,iXi,iEta,iZeta) + U(:, iXi, iEta, EBOTTOM) * spA % bz(iZeta, LEFT) &
                                             * e % geom % normal(IZ,iXi,iEta,EBOTTOM) * e % geom % scal(iXi,iEta,EBOTTOM)
         end do                 ; end do                ; end do

         do iZeta = 0, spA % Nz ; do iEta = 0, spA % Ny ; do iXi = 0, spa % Nx
            faceInt_x(:,iXi,iEta,iZeta) =   faceInt_x(:,iXi,iEta,iZeta) + U(:, iXi, iEta, ETOP) * spA % bz(iZeta, RIGHT) &
                                             * e % geom % normal(IX,iXi,iEta,ETOP) * e % geom % scal(iXi,iEta,ETOP)

            faceInt_y(:,iXi,iEta,iZeta) =   faceInt_y(:,iXi,iEta,iZeta) + U(:, iXi, iEta, ETOP) * spA % bz(iZeta, RIGHT) &
                                             * e % geom % normal(IY,iXi,iEta,ETOP) * e % geom % scal(iXi,iEta,ETOP)

            faceInt_z(:,iXi,iEta,iZeta) =   faceInt_z(:,iXi,iEta,iZeta) + U(:, iXi, iEta, ETOP) * spA % bz(iZeta, RIGHT) &
                                             * e % geom % normal(IZ,iXi,iEta,ETOP) * e % geom % scal(iXi,iEta,ETOP)
         end do                 ; end do                ; end do

      end subroutine VectorWeakIntegrals_StdFace

end module DGWeakIntegrals
