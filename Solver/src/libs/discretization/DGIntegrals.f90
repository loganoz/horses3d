#include "Includes.h"
module DGIntegrals
   use SMConstants
   use ElementClass
   use PhysicsStorage
   use MeshTypes
   use NodalStorageClass, only: NodalStorage

   implicit none

   private

   public ScalarWeakIntegrals_t, VectorWeakIntegrals_t, ScalarStrongIntegrals_t
   public ScalarWeakIntegrals  , VectorWeakIntegrals  , ScalarStrongIntegrals


   type  ScalarWeakIntegrals_t
      contains
         procedure, nopass    :: StdVolumeGreen  => ScalarWeakIntegrals_StdVolumeGreen
         procedure, nopass    :: StdFace => ScalarWeakIntegrals_StdFace
#if defined(NAVIERSTOKES) || defined(INCNS)
         procedure, nopass    :: SplitVolumeDivergence => ScalarWeakIntegrals_SplitVolumeDivergence
#endif
   end type ScalarWeakIntegrals_t

   type  VectorWeakIntegrals_t
      contains
         procedure, nopass    :: StdVolumeGreen  => VectorWeakIntegrals_StdVolumeGreen
         procedure, nopass    :: StdFace => VectorWeakIntegrals_StdFace
   end type VectorWeakIntegrals_t

   type  ScalarStrongIntegrals_t
      contains
         procedure, nopass    :: StdVolumeGreen  => ScalarStrongIntegrals_StdVolumeGreen
#if defined(NAVIERSTOKES) || defined(INCNS)
         procedure, nopass    :: SplitVolumeDivergence => ScalarStrongIntegrals_SplitVolumeDivergence
#endif
   end type ScalarStrongIntegrals_t

   type(ScalarWeakIntegrals_t)   :: ScalarWeakIntegrals
   type(VectorWeakIntegrals_t)   :: VectorWeakIntegrals
   type(ScalarStrongIntegrals_t) :: ScalarStrongIntegrals
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
      function ScalarWeakIntegrals_StdVolumeGreen( e, NEQ, F ) result ( volInt )
         implicit none
         class(Element),      intent(in)  :: e
         integer,             intent(in)  :: NEQ
         real(kind=RP),       intent(in)  :: F     (1:NEQ, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP)                    :: volInt(1:NEQ, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j, k, l
         associate(spAxi   => NodalStorage(e % Nxyz(1)), &
                   spAeta  => NodalStorage(e % Nxyz(2)), &
                   spAzeta => NodalStorage(e % Nxyz(3)) )

         volInt = 0.0_RP

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2)   ; do l = 0, e%Nxyz(1) ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + spAxi % hatD(i,l) * F(:,l,j,k,IX)
         end do             ; end do               ; end do             ; end do

         do k = 0, e%Nxyz(3) ; do l = 0, e%Nxyz(2) ; do j = 0, e%Nxyz(2)   ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + spAeta % hatD(j,l) * F(:,i,l,k,IY)
         end do             ; end do               ; end do             ; end do

         do l = 0, e%Nxyz(3) ; do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2)   ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + spAzeta % hatD(k,l) * F(:,i,j,l,IZ)
         end do             ; end do             ; end do               ; end do

         end associate
      end function ScalarWeakIntegrals_StdVolumeGreen
!
!/////////////////////////////////////////////////////////////////////////////////
!
#if defined(NAVIERSTOKES) || defined(INCNS)
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

         associate(spAxi   => NodalStorage(e % Nxyz(1)), &
                   spAeta  => NodalStorage(e % Nxyz(2)), &
                   spAzeta => NodalStorage(e % Nxyz(3)) )

         volInt = 0.0_RP

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2)   ; do l = 0, e%Nxyz(1) ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + spAxi % sharpD(i,l) * fSharp(:,l,i,j,k) &
                                              + spAxi % hatD(i,l) * Fv(:,l,j,k,IX)
         end do             ; end do               ; end do             ; end do

         do k = 0, e%Nxyz(3) ; do l = 0, e%Nxyz(2) ; do j = 0, e%Nxyz(2)   ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + spAeta % sharpD(j,l) * gSharp(:,l,i,j,k) &
                                              + spAeta % hatD(j,l) * Fv(:,i,l,k,IY)
         end do             ; end do               ; end do             ; end do

         do l = 0, e%Nxyz(3) ; do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2)   ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + spAzeta % sharpD(k,l) * hSharp(:,l,i,j,k) &
                                              + spAzeta % hatD(k,l) * Fv(:,i,j,l,IZ)
         end do             ; end do             ; end do               ; end do

         end associate
      end function ScalarWeakIntegrals_SplitVolumeDivergence
#endif
!
!/////////////////////////////////////////////////////////////////////////////////
!
      pure function ScalarWeakIntegrals_StdFace(e, NEQ, F_FR, F_BK, F_BOT, F_R, F_T, F_L) result(faceInt)
         implicit none
         class(Element),      intent(in)     :: e
         integer,             intent(in)     :: NEQ
         real(kind=RP),       intent(in)     :: F_FR(1:NEQ,0:e % Nxyz(1),0:e % Nxyz(3))
         real(kind=RP),       intent(in)     :: F_BK(1:NEQ,0:e % Nxyz(1),0:e % Nxyz(3))
         real(kind=RP),       intent(in)     :: F_BOT(1:NEQ,0:e % Nxyz(1),0:e % Nxyz(2))
         real(kind=RP),       intent(in)     :: F_R(1:NEQ,0:e % Nxyz(2),0:e % Nxyz(3))
         real(kind=RP),       intent(in)     :: F_T(1:NEQ,0:e % Nxyz(1),0:e % Nxyz(2))
         real(kind=RP),       intent(in)     :: F_L(1:NEQ,0:e % Nxyz(2),0:e % Nxyz(3))
         real(kind=RP)                       :: faceInt(1:NEQ, 0:e%Nxyz(1),0:e%Nxyz(2),0:e%Nxyz(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: iXi, iEta, iZeta, iVar

         associate(spAxi   => NodalStorage(e % Nxyz(1)), &
                   spAeta  => NodalStorage(e % Nxyz(2)), &
                   spAzeta => NodalStorage(e % Nxyz(3)) )
!
!        ----------------
!>       Xi-contributions
!        ----------------
!
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt(:,iXi,iEta,iZeta) = F_L(:, iEta, iZeta) * spAxi % b(iXi, LEFT)
         end do                 ; end do                ; end do

         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt(:,iXi,iEta,iZeta) = faceInt(:,iXi,iEta,iZeta) + F_R(:, iEta, iZeta) * spAxi % b(iXi, RIGHT)
         end do                 ; end do                ; end do
!
!        -----------------
!>       Eta-contributions
!        -----------------
!
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt(:,iXi,iEta,iZeta) = faceInt(:,iXi,iEta,iZeta) + F_FR(:, iXi, iZeta) * spAeta % b(iEta, LEFT)
         end do                 ; end do                ; end do

         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt(:,iXi,iEta,iZeta) = faceInt(:,iXi,iEta,iZeta) + F_BK(:, iXi, iZeta) * spAeta % b(iEta, RIGHT)
         end do                 ; end do                ; end do
!
!        ------------------
!>       Zeta-contributions
!        ------------------
!
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt(:,iXi,iEta,iZeta) = faceInt(:,iXi,iEta,iZeta) + F_BOT(:, iXi, iEta) * spAzeta % b(iZeta, LEFT)
         end do                 ; end do                ; end do

         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt(:,iXi,iEta,iZeta) = faceInt(:,iXi,iEta,iZeta) + F_T(:, iXi, iEta) * spAzeta % b(iZeta, RIGHT)
         end do                 ; end do                ; end do

         end associate
      end function ScalarWeakIntegrals_StdFace
!
!/////////////////////////////////////////////////////////////////////////////
!
!>       Weak integrals with vector test function
!        ----------------------------------------
!/////////////////////////////////////////////////////////////////////////////
!
      subroutine VectorWeakIntegrals_StdVolumeGreen( e, NEQ, U, volInt_x, volInt_y, volInt_z )
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
         integer,             intent(in)  :: NEQ
         real(kind=RP),       intent(in)  :: U        (NEQ,0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),       intent(out) :: volInt_x (NEQ,0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),       intent(out) :: volInt_y (NEQ,0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),       intent(out) :: volInt_z (NEQ,0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i,j,k,l
         real(kind=RP)  :: U_xi(NEQ,0:e % Nxyz(1), 0: e % Nxyz(2), 0: e % Nxyz(3))
         real(kind=RP)  :: U_eta(NEQ,0:e % Nxyz(1), 0: e % Nxyz(2), 0: e % Nxyz(3))
         real(kind=RP)  :: U_zeta(NEQ,0:e % Nxyz(1), 0: e % Nxyz(2), 0: e % Nxyz(3))

         associate(spAxi   => NodalStorage(e % Nxyz(1)), &
                   spAeta  => NodalStorage(e % Nxyz(2)), &
                   spAzeta => NodalStorage(e % Nxyz(3)) )

         volInt_x = 0.0_RP
         volInt_y = 0.0_RP
         volInt_z = 0.0_RP

#ifdef MULTIPHASE
         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2)    ; do l = 0, e%Nxyz(1) ; do i = 0, e%Nxyz(1)
            volInt_x(:,i,j,k) = volInt_x(:,i,j,k) + spAxi % hatD(i,l) * e % geom % jGradXi(IX,l,j,k) * U(:,l,j,k)
            volInt_y(:,i,j,k) = volInt_y(:,i,j,k) + spAxi % hatD(i,l) * e % geom % jGradXi(IY,l,j,k) * U(:,l,j,k)
            volInt_z(:,i,j,k) = volInt_z(:,i,j,k) + spAxi % hatD(i,l) * e % geom % jGradXi(IZ,l,j,k) * U(:,l,j,k)
         end do               ; end do                ; end do             ; end do

         do k = 0, e%Nxyz(3)   ; do l = 0, e%Nxyz(2) ; do j = 0, e%Nxyz(2)    ; do i = 0, e%Nxyz(1)
            volInt_x(:,i,j,k) = volInt_x(:,i,j,k) + spAeta % hatD(j,l) * e % geom % jGradEta(IX,i,l,k) * U(:,i,l,k)
            volInt_y(:,i,j,k) = volInt_y(:,i,j,k) + spAeta % hatD(j,l) * e % geom % jGradEta(IY,i,l,k) * U(:,i,l,k)
            volInt_z(:,i,j,k) = volInt_z(:,i,j,k) + spAeta % hatD(j,l) * e % geom % jGradEta(IZ,i,l,k) * U(:,i,l,k)
         end do               ; end do                ; end do    ; end do

         do l = 0, e%Nxyz(3) ; do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2)    ; do i = 0, e%Nxyz(1)
            volInt_x(:,i,j,k) = volInt_x(:,i,j,k) + spAzeta % hatD(k,l) * e % geom % jGradZeta(IX,i,j,l) * U(:,i,j,l)
            volInt_y(:,i,j,k) = volInt_y(:,i,j,k) + spAzeta % hatD(k,l) * e % geom % jGradZeta(IY,i,j,l) * U(:,i,j,l)
            volInt_z(:,i,j,k) = volInt_z(:,i,j,k) + spAzeta % hatD(k,l) * e % geom % jGradZeta(IZ,i,j,l) * U(:,i,j,l)
         end do             ; end do               ; end do                ; end do
#else
         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2)    ; do l = 0, e%Nxyz(1) ; do i = 0, e%Nxyz(1)
            volInt_x(:,i,j,k) = volInt_x(:,i,j,k) + spAxi % hatD(i,l) * e % geom % jGradXi(IX,i,j,k) * U(:,l,j,k)
            volInt_y(:,i,j,k) = volInt_y(:,i,j,k) + spAxi % hatD(i,l) * e % geom % jGradXi(IY,i,j,k) * U(:,l,j,k)
            volInt_z(:,i,j,k) = volInt_z(:,i,j,k) + spAxi % hatD(i,l) * e % geom % jGradXi(IZ,i,j,k) * U(:,l,j,k)
         end do               ; end do                ; end do             ; end do

         do k = 0, e%Nxyz(3)   ; do l = 0, e%Nxyz(2) ; do j = 0, e%Nxyz(2)    ; do i = 0, e%Nxyz(1)
            volInt_x(:,i,j,k) = volInt_x(:,i,j,k) + spAeta % hatD(j,l) * e % geom % jGradEta(IX,i,j,k) * U(:,i,l,k)
            volInt_y(:,i,j,k) = volInt_y(:,i,j,k) + spAeta % hatD(j,l) * e % geom % jGradEta(IY,i,j,k) * U(:,i,l,k)
            volInt_z(:,i,j,k) = volInt_z(:,i,j,k) + spAeta % hatD(j,l) * e % geom % jGradEta(IZ,i,j,k) * U(:,i,l,k)
         end do               ; end do                ; end do    ; end do

         do l = 0, e%Nxyz(3) ; do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2)    ; do i = 0, e%Nxyz(1)
            volInt_x(:,i,j,k) = volInt_x(:,i,j,k) + spAzeta % hatD(k,l) * e % geom % jGradZeta(IX,i,j,k) * U(:,i,j,l)
            volInt_y(:,i,j,k) = volInt_y(:,i,j,k) + spAzeta % hatD(k,l) * e % geom % jGradZeta(IY,i,j,k) * U(:,i,j,l)
            volInt_z(:,i,j,k) = volInt_z(:,i,j,k) + spAzeta % hatD(k,l) * e % geom % jGradZeta(IZ,i,j,k) * U(:,i,j,l)
         end do             ; end do               ; end do                ; end do

#endif


         end associate
      end subroutine VectorWeakIntegrals_StdVolumeGreen
!
!/////////////////////////////////////////////////////////////////////////////////
!
      subroutine VectorWeakIntegrals_StdFace( e, NEQ, HF, HBK, HBO, HR, HT, HL , &
                                             faceInt_x, faceInt_y, faceInt_z )
         use ElementClass
         use Physics
         use PhysicsStorage
         implicit none
         class(Element), intent(in)  :: e
         integer,        intent(in)  :: NEQ
         real(kind=RP),  intent(in)  :: HF  (NEQ,NDIM, 0:e % Nxyz(1), 0: e % Nxyz(3))
         real(kind=RP),  intent(in)  :: HBK (NEQ,NDIM, 0:e % Nxyz(1), 0: e % Nxyz(3))
         real(kind=RP),  intent(in)  :: HBO (NEQ,NDIM, 0:e % Nxyz(1), 0: e % Nxyz(2))
         real(kind=RP),  intent(in)  :: HR  (NEQ,NDIM, 0:e % Nxyz(2), 0: e % Nxyz(3))
         real(kind=RP),  intent(in)  :: HT  (NEQ,NDIM, 0:e % Nxyz(1), 0: e % Nxyz(2))
         real(kind=RP),  intent(in)  :: HL  (NEQ,NDIM, 0:e % Nxyz(2), 0: e % Nxyz(3))
         real(kind=RP),  intent(out) :: faceInt_x(NEQ, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),  intent(out) :: faceInt_y(NEQ, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),  intent(out) :: faceInt_z(NEQ, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: iVar, iXi, iEta, iZeta

         associate(spAxi   => NodalStorage(e % Nxyz(1)), &
                   spAeta  => NodalStorage(e % Nxyz(2)), &
                   spAzeta => NodalStorage(e % Nxyz(3)) )

!
!        ----------------
!>       Xi-contributions
!        ----------------
!
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt_x(:,iXi,iEta,iZeta) =   HL(:, IX, iEta, iZeta) * spAxi % b(iXi, LEFT)
            faceInt_y(:,iXi,iEta,iZeta) =   HL(:, IY, iEta, iZeta) * spAxi % b(iXi, LEFT)
            faceInt_z(:,iXi,iEta,iZeta) =   HL(:, IZ, iEta, iZeta) * spAxi % b(iXi, LEFT)
         end do                 ; end do                ; end do

         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt_x(:,iXi,iEta,iZeta) =   faceInt_x(:,iXi,iEta,iZeta) + HR(:,IX, iEta, iZeta) * spAxi % b(iXi, RIGHT)
            faceInt_y(:,iXi,iEta,iZeta) =   faceInt_y(:,iXi,iEta,iZeta) + HR(:,IY, iEta, iZeta) * spAxi % b(iXi, RIGHT)
            faceInt_z(:,iXi,iEta,iZeta) =   faceInt_z(:,iXi,iEta,iZeta) + HR(:,IZ, iEta, iZeta) * spAxi % b(iXi, RIGHT)
         end do                 ; end do                ; end do
!
!        -----------------
!>       Eta-contributions
!        -----------------
!
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt_x(:,iXi,iEta,iZeta) =   faceInt_x(:,iXi,iEta,iZeta) + HF(:,IX, iXi, iZeta) * spAeta % b(iEta, LEFT)
            faceInt_y(:,iXi,iEta,iZeta) =   faceInt_y(:,iXi,iEta,iZeta) + HF(:,IY, iXi, iZeta) * spAeta % b(iEta, LEFT)
            faceInt_z(:,iXi,iEta,iZeta) =   faceInt_z(:,iXi,iEta,iZeta) + HF(:,IZ, iXi, iZeta) * spAeta % b(iEta, LEFT)
         end do                 ; end do                ; end do

         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt_x(:,iXi,iEta,iZeta) =   faceInt_x(:,iXi,iEta,iZeta) + HBK(:,IX, iXi, iZeta) * spAeta % b(iEta, RIGHT)
            faceInt_y(:,iXi,iEta,iZeta) =   faceInt_y(:,iXi,iEta,iZeta) + HBK(:,IY, iXi, iZeta) * spAeta % b(iEta, RIGHT)
            faceInt_z(:,iXi,iEta,iZeta) =   faceInt_z(:,iXi,iEta,iZeta) + HBK(:,IZ, iXi, iZeta) * spAeta % b(iEta, RIGHT)
         end do                 ; end do                ; end do
!
!        ------------------
!>       Zeta-contributions
!        ------------------
!
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt_x(:,iXi,iEta,iZeta) =   faceInt_x(:,iXi,iEta,iZeta) + HBO(:, IX, iXi, iEta) * spAzeta % b(iZeta, LEFT)
            faceInt_y(:,iXi,iEta,iZeta) =   faceInt_y(:,iXi,iEta,iZeta) + HBO(:, IY, iXi, iEta) * spAzeta % b(iZeta, LEFT)
            faceInt_z(:,iXi,iEta,iZeta) =   faceInt_z(:,iXi,iEta,iZeta) + HBO(:, IZ, iXi, iEta) * spAzeta % b(iZeta, LEFT)
         end do                 ; end do                ; end do

         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            faceInt_x(:,iXi,iEta,iZeta) =   faceInt_x(:,iXi,iEta,iZeta) + HT(:, IX, iXi, iEta) * spAzeta % b(iZeta, RIGHT)
            faceInt_y(:,iXi,iEta,iZeta) =   faceInt_y(:,iXi,iEta,iZeta) + HT(:, IY, iXi, iEta) * spAzeta % b(iZeta, RIGHT)
            faceInt_z(:,iXi,iEta,iZeta) =   faceInt_z(:,iXi,iEta,iZeta) + HT(:, IZ, iXi, iEta) * spAzeta % b(iZeta, RIGHT)
         end do                 ; end do                ; end do

         end associate
      end subroutine VectorWeakIntegrals_StdFace
!
!/////////////////////////////////////////////////////////////////////////////
!
!>       Strong integrals with scalar test function
!        ------------------------------------------
!/////////////////////////////////////////////////////////////////////////////
!
      function ScalarStrongIntegrals_StdVolumeGreen( e, NEQ, F ) result ( volInt )
         implicit none
         class(Element),      intent(in)  :: e
         integer,             intent(in)  :: NEQ
         real(kind=RP),       intent(in)  :: F     (1:NEQ, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP)                    :: volInt(1:NEQ, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j, k, l

         associate(spAxi   => NodalStorage(e % Nxyz(1)), &
                   spAeta  => NodalStorage(e % Nxyz(2)), &
                   spAzeta => NodalStorage(e % Nxyz(3)) )
         volInt = 0.0_RP

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2)   ; do l = 0, e%Nxyz(1) ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + spAxi % D(i,l) * F(:,l,j,k,IX)
         end do             ; end do               ; end do             ; end do

         do k = 0, e%Nxyz(3) ; do l = 0, e%Nxyz(2) ; do j = 0, e%Nxyz(2)   ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + spAeta % D(j,l) * F(:,i,l,k,IY)
         end do             ; end do               ; end do             ; end do

         do l = 0, e%Nxyz(3) ; do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2)   ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + spAzeta % D(k,l) * F(:,i,j,l,IZ)
         end do             ; end do             ; end do               ; end do

         end associate
      end function ScalarStrongIntegrals_StdVolumeGreen
!
!/////////////////////////////////////////////////////////////////////////////////
!
#if defined(NAVIERSTOKES) || defined(INCNS)
      function ScalarStrongIntegrals_SplitVolumeDivergence( e, fSharp, gSharp, hSharp, Fv ) result ( volInt )
         implicit none
         class(Element), intent(in)  :: e
         real(kind=RP),  intent(in)  :: fSharp(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),  intent(in)  :: gSharp(1:NCONS, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),  intent(in)  :: hSharp(1:NCONS, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),  intent(in)  :: Fv(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP)               :: volInt(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer :: i, j, k, l


         associate(spAxi   => NodalStorage(e % Nxyz(1)), &
                   spAeta  => NodalStorage(e % Nxyz(2)), &
                   spAzeta => NodalStorage(e % Nxyz(3)) )

         volInt = 0.0_RP

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2)   ; do l = 0, e%Nxyz(1) ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + 2.0_RP * spAxi % D(i,l) * fSharp(:,l,i,j,k) &
                                              + spAxi % hatD(i,l) * Fv(:,l,j,k,IX)
         end do             ; end do               ; end do             ; end do

         do k = 0, e%Nxyz(3) ; do l = 0, e%Nxyz(2) ; do j = 0, e%Nxyz(2)   ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + 2.0_RP * spAeta % D(j,l) * gSharp(:,l,i,j,k) &
                                              + spAeta % hatD(j,l) * Fv(:,i,l,k,IY)
         end do             ; end do               ; end do             ; end do

         do l = 0, e%Nxyz(3) ; do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2)   ; do i = 0, e%Nxyz(1)
            volInt(:,i,j,k) = volInt(:,i,j,k) + 2.0_RP * spAzeta % D(k,l) * hSharp(:,l,i,j,k) &
                                              + spAzeta % hatD(k,l) * Fv(:,i,j,l,IZ)
         end do             ; end do             ; end do               ; end do

         end associate

      end function ScalarStrongIntegrals_SplitVolumeDivergence
#endif

end module DGIntegrals
