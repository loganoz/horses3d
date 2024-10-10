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

   public ScalarWeakIntegrals_StdVolumeGreen, ScalarWeakIntegrals_StdFace
   public VectorWeakIntegrals_StdFace


   type  ScalarWeakIntegrals_t
      contains
         procedure, nopass    :: StdVolumeGreen  => ScalarWeakIntegrals_StdVolumeGreen
         procedure, nopass    :: StdFace => ScalarWeakIntegrals_StdFace
#if defined(NAVIERSTOKES) || defined(INCNS)
         procedure, nopass    :: SplitVolumeDivergence => ScalarWeakIntegrals_SplitVolumeDivergence
         procedure, nopass    :: TelescopicVolumeDivergence => ScalarWeakIntegrals_TelescopicVolumeDivergence
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
         procedure, nopass    :: TelescopicVolumeDivergence => ScalarStrongIntegrals_TelescopicVolumeDivergence
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
      subroutine ScalarWeakIntegrals_StdVolumeGreen( Nxyz, NEQ, F, volInt )
         !$acc routine vector
         implicit none
         integer,             intent(in)  :: Nxyz(3)
         integer,             intent(in)  :: NEQ
         real(kind=RP),       intent(in)  :: F     (1:NEQ, 0:Nxyz(1), 0:Nxyz(2), 0:Nxyz(3), 1:NDIM)
         real(kind=RP),    intent(inout)  :: volInt(1:NEQ, 0:Nxyz(1), 0:Nxyz(2), 0:Nxyz(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer   :: i, j, k, l

         volInt = 0.0_RP

         !$acc loop vector collapse(3)
          do k = 0, Nxyz(3) ; do j = 0, Nxyz(2) ; do i = 0, Nxyz(1)  
            !$acc loop seq 
            do l = 0, Nxyz(1)
               volInt(:,i,j,k) = volInt(:,i,j,k) +  NodalStorage(Nxyz(1)) % hatD(i,l) * F(:,l,j,k,IX) &
                                                 +  NodalStorage(Nxyz(2)) % hatD(j,l) * F(:,i,l,k,IY) &
                                                 +  NodalStorage(Nxyz(3)) % hatD(k,l) * F(:,i,j,l,IZ)
            end do             
         end do               ; end do             ; end do

      end subroutine ScalarWeakIntegrals_StdVolumeGreen
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
!
!/////////////////////////////////////////////////////////////////////////////////
!
      function ScalarWeakIntegrals_TelescopicVolumeDivergence(e, Fx, Fy, Fz, Fv) result(volInt)
!
!        ---------
!        Interface
!        ---------
         implicit none
         class(Element), intent(in) :: e
         real(RP),       intent(in) :: Fx    (1:NCONS, 0:e%Nxyz(1)+1, 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(RP),       intent(in) :: Fy    (1:NCONS, 0:e%Nxyz(2)+1, 0:e%Nxyz(1), 0:e%Nxyz(3))
         real(RP),       intent(in) :: Fz    (1:NCONS, 0:e%Nxyz(3)+1, 0:e%Nxyz(1), 0:e%Nxyz(2))
         real(RP),       intent(in) :: Fv    (1:NCONS, 0:e%Nxyz(1),   0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM)
         real(RP)                   :: volInt(1:NCONS, 0:e%Nxyz(1),   0:e%Nxyz(2), 0:e%Nxyz(3))
!
!        ---------------
!        Local variables
!        ---------------
         integer :: i, j, k, l


         associate(Nx => e % Nxyz(1), &
                   Ny => e % Nxyz(2), &
                   Nz => e % Nxyz(3)  )
         associate(spAxi   => NodalStorage(Nx), &
                   spAeta  => NodalStorage(Ny), &
                   spAzeta => NodalStorage(Nz)  )

         volInt = 0.0_RP
!
!        Xi
!        --
         do k = 0, Nz ; do j = 0, Ny ; do i = 0, Nx
            volInt(:,i,j,k) = volInt(:,i,j,k) + (Fx(:,i+1,j,k)-Fx(:,i,j,k)) / spAxi % w(i)
         end do       ; end do       ; end do

         do k = 0, Nz ; do j = 0, Ny
            volInt(:,0,j,k)  = volInt(:,0,j,k)  + Fx(:,0,j,k)    / spAxi % w(0)
            volInt(:,Nx,j,k) = volInt(:,Nx,j,k) - Fx(:,Nx+1,j,k) / spAxi % w(Nx)
         end do       ; end do

         do k = 0, Nz ; do j = 0, Ny ; do l = 0, Nx ; do i = 0, Nx
            volInt(:,i,j,k) = volInt(:,i,j,k) + spAxi % hatD(i,l) * Fv(:,l,j,k,IX)
         end do       ; end do       ; end do       ; end do
!
!        Eta
!        ---
         do k = 0, Nz ; do i = 0, Nx ; do j = 0, Ny
            volInt(:,i,j,k) = volInt(:,i,j,k) + (Fy(:,j+1,i,k)-Fy(:,j,i,k)) / spAeta % w(j)
         end do       ; end do       ; end do

         do k = 0, Nz ; do i = 0, Nx
            volInt(:,i,0,k)  = volInt(:,i,0,k)  + Fy(:,0,i,k)    / spAeta % w(0)
            volInt(:,i,Ny,k) = volInt(:,i,Ny,k) - Fy(:,Ny+1,i,k) / spAeta % w(Ny)
         end do       ; end do

         do k = 0, Nz ; do l = 0, Ny ; do j = 0, Ny ; do i = 0, Nx
            volInt(:,i,j,k) = volInt(:,i,j,k) + spAeta % hatD(j,l) * Fv(:,i,l,k,IY)
         end do       ; end do       ; end do       ; end do
!
!        Zeta
!        ----
         do j = 0, Ny ; do i = 0, Nx ; do k = 0, Nz
            volInt(:,i,j,k) = volInt(:,i,j,k) + (Fz(:,k+1,i,j)-Fz(:,k,i,j)) / spAzeta % w(k)
         end do       ; end do       ; end do

         do j = 0, Ny ; do i = 0, Nx
            volInt(:,i,j,0)  = volInt(:,i,j,0)  + Fz(:,0,i,j)    / spAzeta % w(0)
            volInt(:,i,j,Nz) = volInt(:,i,j,Nz) - Fz(:,Nz+1,i,j) / spAzeta % w(Nz)
         end do       ; end do

         do l = 0, Nz ; do k = 0, Nz ; do j = 0, Ny ; do i = 0, Nx
            volInt(:,i,j,k) = volInt(:,i,j,k) + spAzeta % hatD(k,l) * Fv(:,i,j,l,IZ)
         end do       ; end do       ; end do       ; end do

         end associate
         end associate

      end function ScalarWeakIntegrals_TelescopicVolumeDivergence
#endif
!
!///////////////////////////////////////////////////////////////
!
      subroutine ScalarWeakIntegrals_StdFace(NEQ, Nxyz, F_FR, F_BK, F_BOT, F_R, F_T, F_L, intFace)
         !$acc routine vector
         implicit none
         integer,             intent(in)     :: NEQ
         integer,             intent(in)     :: Nxyz(3)
         real(kind=RP),       intent(in)     :: F_FR(1:NEQ,0:Nxyz(1),0:Nxyz(3))
         real(kind=RP),       intent(in)     :: F_BK(1:NEQ,0:Nxyz(1),0:Nxyz(3))
         real(kind=RP),       intent(in)     :: F_BOT(1:NEQ,0:Nxyz(1),0:Nxyz(2))
         real(kind=RP),       intent(in)     :: F_R(1:NEQ,0:Nxyz(2),0:Nxyz(3))
         real(kind=RP),       intent(in)     :: F_T(1:NEQ,0:Nxyz(1),0:Nxyz(2))
         real(kind=RP),       intent(in)     :: F_L(1:NEQ,0:Nxyz(2),0:Nxyz(3))
         real(kind=RP),      intent(inout)   :: intFace(1:NEQ, 0:Nxyz(1),0:Nxyz(2),0:Nxyz(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: iXi, iEta, iZeta,eq
  
         !$acc loop vector collapse(3)
         do iZeta = 0, Nxyz(3) 
            do iEta = 0, Nxyz(2) 
               do iXi = 0, Nxyz(1)
                  !$acc loop seq
                  do eq = 1, NCONS
                     intFace(eq,iXi,iEta,iZeta) = intFace(eq,iXi,iEta,iZeta) - ( &
                                                + F_L(eq, iEta, iZeta) * NodalStorage(Nxyz(1)) % b(iXi, LEFT)    &
                                                + F_R(eq, iEta, iZeta) * NodalStorage(Nxyz(1)) % b(iXi, RIGHT)   &
                                                + F_FR(eq, iXi, iZeta) * NodalStorage(Nxyz(2)) % b(iEta, LEFT)   &
                                                + F_BK(eq, iXi, iZeta) * NodalStorage(Nxyz(2)) % b(iEta, RIGHT)  &
                                                + F_BOT(eq, iXi, iEta) * NodalStorage(Nxyz(3)) % b(iZeta, LEFT)  &
                                                + F_T(eq, iXi, iEta)   * NodalStorage(Nxyz(3)) % b(iZeta, RIGHT) ) 
                     enddo
               end do                 
            end do                
         end do

      end subroutine ScalarWeakIntegrals_StdFace
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
         !$acc routine vector 
         use ElementClass
         use Physics
         use PhysicsStorage
         implicit none
         type(Element), intent(in)  :: e
         integer,        intent(in)  :: NEQ
         real(kind=RP),  intent(in)  :: HF  (NEQ,NDIM, 0:e % Nxyz(1), 0: e % Nxyz(3))
         real(kind=RP),  intent(in)  :: HBK (NEQ,NDIM, 0:e % Nxyz(1), 0: e % Nxyz(3))
         real(kind=RP),  intent(in)  :: HBO (NEQ,NDIM, 0:e % Nxyz(1), 0: e % Nxyz(2))
         real(kind=RP),  intent(in)  :: HR  (NEQ,NDIM, 0:e % Nxyz(2), 0: e % Nxyz(3))
         real(kind=RP),  intent(in)  :: HT  (NEQ,NDIM, 0:e % Nxyz(1), 0: e % Nxyz(2))
         real(kind=RP),  intent(in)  :: HL  (NEQ,NDIM, 0:e % Nxyz(2), 0: e % Nxyz(3))
         real(kind=RP),  intent(inout) :: faceInt_x(NEQ, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),  intent(inout) :: faceInt_y(NEQ, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP),  intent(inout) :: faceInt_z(NEQ, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: iVar, iXi, iEta, iZeta, eq
!
!        ----------------
!>       Xi-contributions
!        ----------------
!
         !$acc loop vector collapse(3)
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            !$acc loop seq
            do eq = 1, NCONS
               faceInt_x(eq,iXi,iEta,iZeta) =   faceInt_x(eq,iXi,iEta,iZeta) + HL(eq, IX, iEta, iZeta) * NodalStorage(e % Nxyz(1)) % b(iXi, LEFT)
               faceInt_y(eq,iXi,iEta,iZeta) =   faceInt_y(eq,iXi,iEta,iZeta) + HL(eq, IY, iEta, iZeta) * NodalStorage(e % Nxyz(1)) % b(iXi, LEFT)
               faceInt_z(eq,iXi,iEta,iZeta) =   faceInt_z(eq,iXi,iEta,iZeta) + HL(eq, IZ, iEta, iZeta) * NodalStorage(e % Nxyz(1)) % b(iXi, LEFT)
            enddo
         end do                 ; end do                ; end do
         
         !$acc loop vector collapse(3)
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            !$acc loop seq
            do eq = 1, NCONS
               faceInt_x(eq,iXi,iEta,iZeta) =   faceInt_x(eq,iXi,iEta,iZeta) + HR(eq,IX, iEta, iZeta) * NodalStorage(e % Nxyz(1)) % b(iXi, RIGHT)
               faceInt_y(eq,iXi,iEta,iZeta) =   faceInt_y(eq,iXi,iEta,iZeta) + HR(eq,IY, iEta, iZeta) * NodalStorage(e % Nxyz(1)) % b(iXi, RIGHT)
               faceInt_z(eq,iXi,iEta,iZeta) =   faceInt_z(eq,iXi,iEta,iZeta) + HR(eq,IZ, iEta, iZeta) * NodalStorage(e % Nxyz(1)) % b(iXi, RIGHT)
            enddo
         end do                 ; end do                ; end do
!
!        -----------------
!>       Eta-contributions
!        -----------------
!
         !$acc loop vector collapse(3)
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            !$acc loop seq
            do eq = 1, NCONS
               faceInt_x(eq,iXi,iEta,iZeta) =   faceInt_x(eq,iXi,iEta,iZeta) + HF(eq,IX, iXi, iZeta) * NodalStorage(e % Nxyz(2)) % b(iEta, LEFT)
               faceInt_y(eq,iXi,iEta,iZeta) =   faceInt_y(eq,iXi,iEta,iZeta) + HF(eq,IY, iXi, iZeta) * NodalStorage(e % Nxyz(2)) % b(iEta, LEFT)
               faceInt_z(eq,iXi,iEta,iZeta) =   faceInt_z(eq,iXi,iEta,iZeta) + HF(eq,IZ, iXi, iZeta) * NodalStorage(e % Nxyz(2)) % b(iEta, LEFT)
            enddo
         end do                 ; end do                ; end do

         !$acc loop vector collapse(3)
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            !$acc loop seq
            do eq = 1, NCONS
               faceInt_x(eq,iXi,iEta,iZeta) =   faceInt_x(eq,iXi,iEta,iZeta) + HBK(eq,IX, iXi, iZeta) * NodalStorage(e % Nxyz(2)) % b(iEta, RIGHT)
               faceInt_y(eq,iXi,iEta,iZeta) =   faceInt_y(eq,iXi,iEta,iZeta) + HBK(eq,IY, iXi, iZeta) * NodalStorage(e % Nxyz(2)) % b(iEta, RIGHT)
               faceInt_z(eq,iXi,iEta,iZeta) =   faceInt_z(eq,iXi,iEta,iZeta) + HBK(eq,IZ, iXi, iZeta) * NodalStorage(e % Nxyz(2)) % b(iEta, RIGHT)
            enddo
         end do                 ; end do                ; end do
!
!        ------------------
!>       Zeta-contributions
!        ------------------
!
         !$acc loop vector collapse(3)
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            !$acc loop seq
            do eq = 1, NCONS
               faceInt_x(eq,iXi,iEta,iZeta) =   faceInt_x(eq,iXi,iEta,iZeta) + HBO(eq, IX, iXi, iEta) * NodalStorage(e % Nxyz(3)) % b(iZeta, LEFT)
               faceInt_y(eq,iXi,iEta,iZeta) =   faceInt_y(eq,iXi,iEta,iZeta) + HBO(eq, IY, iXi, iEta) * NodalStorage(e % Nxyz(3)) % b(iZeta, LEFT)
               faceInt_z(eq,iXi,iEta,iZeta) =   faceInt_z(eq,iXi,iEta,iZeta) + HBO(eq, IZ, iXi, iEta) * NodalStorage(e % Nxyz(3)) % b(iZeta, LEFT)
            enddo
         end do                 ; end do                ; end do
         
         !$acc loop vector collapse(3)
         do iZeta = 0, e%Nxyz(3) ; do iEta = 0, e%Nxyz(2) ; do iXi = 0, e%Nxyz(1)
            !$acc loop seq
            do eq = 1, NCONS
               faceInt_x(eq,iXi,iEta,iZeta) =   faceInt_x(eq,iXi,iEta,iZeta) + HT(eq, IX, iXi, iEta) * NodalStorage(e % Nxyz(3)) % b(iZeta, RIGHT)
               faceInt_y(eq,iXi,iEta,iZeta) =   faceInt_y(eq,iXi,iEta,iZeta) + HT(eq, IY, iXi, iEta) * NodalStorage(e % Nxyz(3)) % b(iZeta, RIGHT)
               faceInt_z(eq,iXi,iEta,iZeta) =   faceInt_z(eq,iXi,iEta,iZeta) + HT(eq, IZ, iXi, iEta) * NodalStorage(e % Nxyz(3)) % b(iZeta, RIGHT)
            enddo
         end do                 ; end do                ; end do

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
!
!/////////////////////////////////////////////////////////////////////////////////
!
      function ScalarStrongIntegrals_TelescopicVolumeDivergence(e, Fx, Fy, Fz, Fv) result(volInt)
!
!        ---------
!        Interface
!        ---------
         implicit none
         class(Element), intent(in) :: e
         real(RP),       intent(in) :: Fx    (1:NCONS, 0:e%Nxyz(1)+1, 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(RP),       intent(in) :: Fy    (1:NCONS, 0:e%Nxyz(2)+1, 0:e%Nxyz(1), 0:e%Nxyz(3))
         real(RP),       intent(in) :: Fz    (1:NCONS, 0:e%Nxyz(3)+1, 0:e%Nxyz(1), 0:e%Nxyz(2))
         real(RP),       intent(in) :: Fv    (1:NCONS, 0:e%Nxyz(1),   0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM)
         real(RP)                   :: volInt(1:NCONS, 0:e%Nxyz(1),   0:e%Nxyz(2), 0:e%Nxyz(3))
!
!        ---------------
!        Local variables
!        ---------------
         integer :: i, j, k, l


         associate(Nx => e % Nxyz(1), &
                   Ny => e % Nxyz(2), &
                   Nz => e % Nxyz(3)  )
         associate(spAxi   => NodalStorage(Nx), &
                   spAeta  => NodalStorage(Ny), &
                   spAzeta => NodalStorage(Nz)  )

         volInt = 0.0_RP
!
!        Xi
!        --
         do k = 0, Nz ; do j = 0, Ny ; do i = 0, Nx
            volInt(:,i,j,k) = volInt(:,i,j,k) + (Fx(:,i+1,j,k)-Fx(:,i,j,k)) / spAxi % w(i)
         end do       ; end do       ; end do

         do k = 0, Nz ; do j = 0, Ny ; do l = 0, Nx ; do i = 0, Nx
            volInt(:,i,j,k) = volInt(:,i,j,k) + spAxi % hatD(i,l) * Fv(:,l,j,k,IX)
         end do       ; end do       ; end do       ; end do
!
!        Eta
!        ---
         do k = 0, Nz ; do i = 0, Nx ; do j = 0, Ny
            volInt(:,i,j,k) = volInt(:,i,j,k) + (Fy(:,j+1,i,k)-Fy(:,j,i,k)) / spAeta % w(j)
         end do       ; end do       ; end do

         do k = 0, Nz ; do l = 0, Ny ; do j = 0, Ny ; do i = 0, Nx
            volInt(:,i,j,k) = volInt(:,i,j,k) + spAeta % hatD(j,l) * Fv(:,i,l,k,IY)
         end do       ; end do       ; end do       ; end do
!
!        Zeta
!        ----
         do j = 0, Ny ; do i = 0, Nx ; do k = 0, Nz
            volInt(:,i,j,k) = volInt(:,i,j,k) + (Fz(:,k+1,i,j)-Fz(:,k,i,j)) / spAzeta % w(k)
         end do       ; end do       ; end do

         do l = 0, Nz ; do k = 0, Nz ; do j = 0, Ny ; do i = 0, Nx
            volInt(:,i,j,k) = volInt(:,i,j,k) + spAzeta % hatD(k,l) * Fv(:,i,j,l,IZ)
         end do       ; end do       ; end do       ; end do

         end associate
         end associate

      end function ScalarStrongIntegrals_TelescopicVolumeDivergence
#endif

end module DGIntegrals