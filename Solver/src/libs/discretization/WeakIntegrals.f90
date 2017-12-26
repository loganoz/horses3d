module DGWeakIntegrals
   use SMConstants
   use ElementClass
   use PhysicsStorage
   use PhysicsStorage, only: N_EQN, N_GRAD_EQN, IX,IY,IZ
   use MeshTypes
   use VariableConversion, only: gradientValuesForQ
   implicit none


   private

   public ScalarWeakIntegrals_t, VectorWeakIntegrals_t
   public ScalarWeakIntegrals  , VectorWeakIntegrals  
   
   type  ScalarWeakIntegrals_t
      contains
         procedure, nopass    :: StdVolumeGreen  => ScalarWeakIntegrals_StdVolumeGreen
         procedure, nopass    :: SplitVolumeDivergence => ScalarWeakIntegrals_SplitVolumeDivergence
         procedure, nopass    :: StdFace => ScalarWeakIntegrals_StdFace
   end type ScalarWeakIntegrals_t

   type  VectorWeakIntegrals_t
      contains
         procedure, nopass    :: StdVolumeGreen  => VectorWeakIntegrals_StdVolumeGreen
         procedure, nopass    :: StdFace => VectorWeakIntegrals_StdFace
   end type VectorWeakIntegrals_t


   type(ScalarWeakIntegrals_t)   :: ScalarWeakIntegrals
   type(VectorWeakIntegrals_t)   :: VectorWeakIntegrals
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
      function ScalarWeakIntegrals_StdVolumeGreen( e, F ) result ( volInt )
         implicit none
         class(Element),      intent(in)  :: e
         real(kind=RP),       intent(in)  :: F     (1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP)                    :: volInt(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j, k, l

         volInt = 0.0_RP

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2)   ; do l = 0, e%Nxyz(1) ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + e % spAxi % hatD(i,l) * F(:,l,j,k,IX)
         end do             ; end do               ; end do             ; end do

         do k = 0, e%Nxyz(3) ; do l = 0, e%Nxyz(2) ; do j = 0, e%Nxyz(2)   ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + e % spAeta % hatD(j,l) * F(:,i,l,k,IY)
         end do             ; end do               ; end do             ; end do

         do l = 0, e%Nxyz(3) ; do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2)   ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + e % spAzeta % hatD(k,l) * F(:,i,j,l,IZ)
         end do             ; end do             ; end do               ; end do

      end function ScalarWeakIntegrals_StdVolumeGreen

      function ScalarWeakIntegrals_SplitVolumeDivergence( e, fSharp, gSharp, hSharp, Fv ) result ( volInt )
         implicit none
         class(Element),      intent(in)  :: e
         real(kind=RP),       intent(in)  :: fSharp(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),       intent(in)  :: gSharp(1:NCONS, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),       intent(in)  :: hSharp(1:NCONS, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),       intent(in)  :: Fv(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP)                    :: volInt(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j, k, l

         volInt = 0.0_RP

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2)   ; do l = 0, e%Nxyz(1) ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + e % spAxi % sharpD(i,l) * fSharp(:,l,i,j,k) &
                                              + e % spAxi % hatD(i,l) * Fv(:,l,j,k,IX)
         end do             ; end do               ; end do             ; end do

         do k = 0, e%Nxyz(3) ; do l = 0, e%Nxyz(2) ; do j = 0, e%Nxyz(2)   ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + e % spAeta % sharpD(j,l) * gSharp(:,l,i,j,k) &
                                              + e % spAeta % hatD(j,l) * Fv(:,i,l,k,IY)
         end do             ; end do               ; end do             ; end do

         do l = 0, e%Nxyz(3) ; do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2)   ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + e % spAzeta % sharpD(k,l) * hSharp(:,l,i,j,k) &
                                              + e % spAzeta % hatD(k,l) * Fv(:,i,j,l,IZ)
         end do             ; end do             ; end do               ; end do

      end function ScalarWeakIntegrals_SplitVolumeDivergence

!
!///////////////////////////////////////////////////////////////
!
      pure function ScalarWeakIntegrals_StdFace(e, F_FR, F_BK, F_BOT, F_R, F_T, F_L) result(faceInt)
         implicit none
         class(Element),      intent(in)     :: e
         real(kind=RP),       intent(in)     :: F_FR(1:NCONS,0:e % Nxyz(1),0:e % Nxyz(3))
         real(kind=RP),       intent(in)     :: F_BK(1:NCONS,0:e % Nxyz(1),0:e % Nxyz(3))
         real(kind=RP),       intent(in)     :: F_BOT(1:NCONS,0:e % Nxyz(1),0:e % Nxyz(2))
         real(kind=RP),       intent(in)     :: F_R(1:NCONS,0:e % Nxyz(2),0:e % Nxyz(3))
         real(kind=RP),       intent(in)     :: F_T(1:NCONS,0:e % Nxyz(1),0:e % Nxyz(2))
         real(kind=RP),       intent(in)     :: F_L(1:NCONS,0:e % Nxyz(2),0:e % Nxyz(3))
         real(kind=RP)                       :: faceInt(1:NCONS, 0:e%Nxyz(1),0:e%Nxyz(2),0:e%Nxyz(3))
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
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt(:,iXi,iEta,iZeta) = F_L(:, iEta, iZeta) * e % spAxi % b(iXi, LEFT)
         end do                 ; end do                ; end do
      
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt(:,iXi,iEta,iZeta) = faceInt(:,iXi,iEta,iZeta) + F_R(:, iEta, iZeta) * e % spAxi % b(iXi, RIGHT)
         end do                 ; end do                ; end do
!
!        -----------------
!>       Eta-contributions
!        -----------------
!
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt(:,iXi,iEta,iZeta) = faceInt(:,iXi,iEta,iZeta) + F_FR(:, iXi, iZeta) * e % spAeta % b(iEta, LEFT)
         end do                 ; end do                ; end do

         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt(:,iXi,iEta,iZeta) = faceInt(:,iXi,iEta,iZeta) + F_BK(:, iXi, iZeta) * e % spAeta % b(iEta, RIGHT)
         end do                 ; end do                ; end do
!
!        ------------------
!>       Zeta-contributions
!        ------------------
!
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt(:,iXi,iEta,iZeta) = faceInt(:,iXi,iEta,iZeta) + F_BOT(:, iXi, iEta) * e % spAzeta % b(iZeta, LEFT)
         end do                 ; end do                ; end do

         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt(:,iXi,iEta,iZeta) = faceInt(:,iXi,iEta,iZeta) + F_T(:, iXi, iEta) * e % spAzeta % b(iZeta, RIGHT)
         end do                 ; end do                ; end do

      end function ScalarWeakIntegrals_StdFace
!
!/////////////////////////////////////////////////////////////////////////////
!
!>       Weak integrals with vector test function
!        ----------------------------------------
!/////////////////////////////////////////////////////////////////////////////
!
      subroutine VectorWeakIntegrals_StdVolumeGreen( e, U, volInt_x, volInt_y, volInt_z )
!
!        ***********************************************************************************
!              This integrals compute:
!
!                 volInt_d(i,j,k) =   (Ja^1_d)_{ijk}\sum_{m=0}^N \hat{D}_{im}U_{mjk}
!                                   + (Ja^2_d)_{ijk}\sum_{m=0}^N \hat{D}_{jm}U_{imk}
!                                   + (Ja^3_d)_{ijk}\sum_{m=0}^N \hat{D}_{km}U_{ijk}
!
!        ***********************************************************************************
!
         use ElementClass
         use Physics
         use PhysicsStorage
         implicit none
         class(Element),      intent(in)  :: e
         real(kind=RP),       intent(in)  :: U        (N_GRAD_EQN,0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),       intent(out) :: volInt_x (N_GRAD_EQN,0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),       intent(out) :: volInt_y (N_GRAD_EQN,0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),       intent(out) :: volInt_z (N_GRAD_EQN,0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i,j,k,l
         real(kind=RP)  :: U_xi(N_GRAD_EQN,0:e % Nxyz(1), 0: e % Nxyz(2), 0: e % Nxyz(3))
         real(kind=RP)  :: U_eta(N_GRAD_EQN,0:e % Nxyz(1), 0: e % Nxyz(2), 0: e % Nxyz(3))
         real(kind=RP)  :: U_zeta(N_GRAD_EQN,0:e % Nxyz(1), 0: e % Nxyz(2), 0: e % Nxyz(3))

         volInt_x = 0.0_RP
         volInt_y = 0.0_RP
         volInt_z = 0.0_RP

         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2)    ; do l = 0, e%Nxyz(1) ; do i = 0, e%Nxyz(1)
            volInt_x(:,i,j,k) = volInt_x(:,i,j,k) + e % spAxi % hatD(i,l) * e % geom % jGradXi(IX,l,j,k) * U(:,l,j,k)
            volInt_y(:,i,j,k) = volInt_y(:,i,j,k) + e % spAxi % hatD(i,l) * e % geom % jGradXi(IY,l,j,k) * U(:,l,j,k)
            volInt_z(:,i,j,k) = volInt_z(:,i,j,k) + e % spAxi % hatD(i,l) * e % geom % jGradXi(IZ,l,j,k) * U(:,l,j,k)
         end do               ; end do                ; end do             ; end do

         do k = 0, e%Nxyz(3)   ; do l = 0, e%Nxyz(2) ; do j = 0, e%Nxyz(2)    ; do i = 0, e%Nxyz(1)
            volInt_x(:,i,j,k) = volInt_x(:,i,j,k) + e % spAeta % hatD(j,l) * e % geom % jGradEta(IX,i,l,k) * U(:,i,l,k)
            volInt_y(:,i,j,k) = volInt_y(:,i,j,k) + e % spAeta % hatD(j,l) * e % geom % jGradEta(IY,i,l,k) * U(:,i,l,k)
            volInt_z(:,i,j,k) = volInt_z(:,i,j,k) + e % spAeta % hatD(j,l) * e % geom % jGradEta(IZ,i,l,k) * U(:,i,l,k)
         end do               ; end do                ; end do    ; end do
 
         do l = 0, e%Nxyz(3) ; do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2)    ; do i = 0, e%Nxyz(1)
            volInt_x(:,i,j,k) = volInt_x(:,i,j,k) + e % spAzeta % hatD(k,l) * e % geom % jGradZeta(IX,i,j,l) * U(:,i,j,l)
            volInt_y(:,i,j,k) = volInt_y(:,i,j,k) + e % spAzeta % hatD(k,l) * e % geom % jGradZeta(IY,i,j,l) * U(:,i,j,l)
            volInt_z(:,i,j,k) = volInt_z(:,i,j,k) + e % spAzeta % hatD(k,l) * e % geom % jGradZeta(IZ,i,j,l) * U(:,i,j,l)
         end do             ; end do               ; end do                ; end do
   
      end subroutine VectorWeakIntegrals_StdVolumeGreen
!
!/////////////////////////////////////////////////////////////////////////////////
!
      subroutine VectorWeakIntegrals_StdFace( e, HF, HBK, HBO, HR, HT, HL , &
                                             faceInt_x, faceInt_y, faceInt_z )
         use ElementClass
         use Physics
         use PhysicsStorage
         implicit none
         class(Element), intent(in)  :: e
         real(kind=RP),  intent(in)  :: HF  (N_GRAD_EQN,NDIM, 0:e % Nxyz(1), 0: e % Nxyz(3))
         real(kind=RP),  intent(in)  :: HBK (N_GRAD_EQN,NDIM, 0:e % Nxyz(1), 0: e % Nxyz(3))
         real(kind=RP),  intent(in)  :: HBO (N_GRAD_EQN,NDIM, 0:e % Nxyz(1), 0: e % Nxyz(2))
         real(kind=RP),  intent(in)  :: HR  (N_GRAD_EQN,NDIM, 0:e % Nxyz(2), 0: e % Nxyz(3))
         real(kind=RP),  intent(in)  :: HT  (N_GRAD_EQN,NDIM, 0:e % Nxyz(1), 0: e % Nxyz(2))
         real(kind=RP),  intent(in)  :: HL  (N_GRAD_EQN,NDIM, 0:e % Nxyz(2), 0: e % Nxyz(3))
         real(kind=RP),  intent(out) :: faceInt_x(N_GRAD_EQN, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3)) 
         real(kind=RP),  intent(out) :: faceInt_y(N_GRAD_EQN, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3)) 
         real(kind=RP),  intent(out) :: faceInt_z(N_GRAD_EQN, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: iVar, iXi, iEta, iZeta
!
!        ----------------
!>       Xi-contributions
!        ----------------
!
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt_x(:,iXi,iEta,iZeta) =   HL(:, IX, iEta, iZeta) * e % spAxi % b(iXi, LEFT)
            faceInt_y(:,iXi,iEta,iZeta) =   HL(:, IY, iEta, iZeta) * e % spAxi % b(iXi, LEFT)
            faceInt_z(:,iXi,iEta,iZeta) =   HL(:, IZ, iEta, iZeta) * e % spAxi % b(iXi, LEFT) 
         end do                 ; end do                ; end do

         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt_x(:,iXi,iEta,iZeta) =   faceInt_x(:,iXi,iEta,iZeta) + HR(:,IX, iEta, iZeta) * e % spAxi % b(iXi, RIGHT)
            faceInt_y(:,iXi,iEta,iZeta) =   faceInt_y(:,iXi,iEta,iZeta) + HR(:,IY, iEta, iZeta) * e % spAxi % b(iXi, RIGHT) 
            faceInt_z(:,iXi,iEta,iZeta) =   faceInt_z(:,iXi,iEta,iZeta) + HR(:,IZ, iEta, iZeta) * e % spAxi % b(iXi, RIGHT) 
         end do                 ; end do                ; end do
!
!        -----------------
!>       Eta-contributions
!        -----------------
!
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt_x(:,iXi,iEta,iZeta) =   faceInt_x(:,iXi,iEta,iZeta) + HF(:,IX, iXi, iZeta) * e % spAeta % b(iEta, LEFT) 
            faceInt_y(:,iXi,iEta,iZeta) =   faceInt_y(:,iXi,iEta,iZeta) + HF(:,IY, iXi, iZeta) * e % spAeta % b(iEta, LEFT) 
            faceInt_z(:,iXi,iEta,iZeta) =   faceInt_z(:,iXi,iEta,iZeta) + HF(:,IZ, iXi, iZeta) * e % spAeta % b(iEta, LEFT) 
         end do                 ; end do                ; end do

         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt_x(:,iXi,iEta,iZeta) =   faceInt_x(:,iXi,iEta,iZeta) + HBK(:,IX, iXi, iZeta) * e % spAeta % b(iEta, RIGHT) 
            faceInt_y(:,iXi,iEta,iZeta) =   faceInt_y(:,iXi,iEta,iZeta) + HBK(:,IY, iXi, iZeta) * e % spAeta % b(iEta, RIGHT) 
            faceInt_z(:,iXi,iEta,iZeta) =   faceInt_z(:,iXi,iEta,iZeta) + HBK(:,IZ, iXi, iZeta) * e % spAeta % b(iEta, RIGHT) 
         end do                 ; end do                ; end do
!
!        ------------------
!>       Zeta-contributions
!        ------------------
!
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt_x(:,iXi,iEta,iZeta) =   faceInt_x(:,iXi,iEta,iZeta) + HBO(:, IX, iXi, iEta) * e % spAzeta % b(iZeta, LEFT) 
            faceInt_y(:,iXi,iEta,iZeta) =   faceInt_y(:,iXi,iEta,iZeta) + HBO(:, IY, iXi, iEta) * e % spAzeta % b(iZeta, LEFT) 
            faceInt_z(:,iXi,iEta,iZeta) =   faceInt_z(:,iXi,iEta,iZeta) + HBO(:, IZ, iXi, iEta) * e % spAzeta % b(iZeta, LEFT) 
         end do                 ; end do                ; end do

         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt_x(:,iXi,iEta,iZeta) =   faceInt_x(:,iXi,iEta,iZeta) + HT(:, IX, iXi, iEta) * e % spAzeta % b(iZeta, RIGHT) 
            faceInt_y(:,iXi,iEta,iZeta) =   faceInt_y(:,iXi,iEta,iZeta) + HT(:, IY, iXi, iEta) * e % spAzeta % b(iZeta, RIGHT) 
            faceInt_z(:,iXi,iEta,iZeta) =   faceInt_z(:,iXi,iEta,iZeta) + HT(:, IZ, iXi, iEta) * e % spAzeta % b(iZeta, RIGHT) 
         end do                 ; end do                ; end do

      end subroutine VectorWeakIntegrals_StdFace
end module DGWeakIntegrals
