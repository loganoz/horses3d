!
!//////////////////////////////////////////////////////
!
!   @File:    Physics_CH.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Thu Apr 19 17:24:31 2018
!   @Last revision date: Tue Apr 24 11:13:20 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 60804273321199c0675663c7d4f1c517987552a7
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module Physics_CH 
!
      USE SMConstants
      USE PhysicsStorage_CH
      IMPLICIT NONE

      private
      public  CHDivergenceFlux, AddQuarticDWPDerivative, QuarticDWP
      public  CHDivergenceFlux0D, CHDivergenceFlux3D, PoiseuilleFlow
!
!     ---------
!     Constants
!     ---------
!
      INTEGER, PARAMETER   :: WALL_BC = 1, RADIATION_BC = 2
      INTEGER              :: boundaryCondition(4), bcType

     interface CHDivergenceFlux
      module procedure CHDivergenceFlux0D, CHDivergenceFlux2D, CHDivergenceFlux3D
     end interface CHDivergenceFlux
!
!     ========
      CONTAINS 
!     ========
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
!>        VISCOUS FLUXES
!         --------------
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
      pure subroutine CHDivergenceFlux0D(nEqn, nGradEqn, Q, U_x, U_y, U_z, mu, kappa, F)
         implicit none
         integer,       intent(in)  :: nEqn, nGradEqn
         real(kind=RP), intent(in)  :: Q   (nEqn)
         real(kind=RP), intent(in)  :: U_x (nGradEqn)
         real(kind=RP), intent(in)  :: U_y (nGradEqn)
         real(kind=RP), intent(in)  :: U_z (nGradEqn)
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: kappa
         real(kind=RP), intent(out) :: F(1:nEqn, 1:NDIM)

         F(1,IX) = U_x(1)
         F(1,IY) = U_y(1)
         F(1,IZ) = U_z(1)

      end subroutine CHDivergenceFlux0D

      pure subroutine CHDivergenceFlux2D(nEqn, nGradEqn, N, Q, U_x, U_y, U_z, mu, kappa, F)
         implicit none
         integer,          intent(in)  :: nEqn, nGradEqn
         integer         , intent(in)  :: N(2)
         real(kind=RP),    intent(in)  :: Q  (1:nEqn, 0:N(1), 0:N(2))
         real(kind=RP),    intent(in)  :: U_x(1:nGradEqn, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: U_y(1:nGradEqn, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: U_z(1:nGradEqn, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: mu  (0:N(1), 0:N(2))
         real(kind=RP),    intent(in)  :: kappa(0:N(1), 0:N(2))
         real(kind=RP),    intent(out) :: F   (1:nEqn, 0:N(1), 0:N(2), 1:NDIM)

         F(1,:,:,IX) = U_x(1,:,:)
         F(1,:,:,IY) = U_y(1,:,:)
         F(1,:,:,IZ) = U_z(1,:,:)

      end subroutine CHDivergenceFlux2D

      pure subroutine CHDivergenceFlux3D(nEqn, nGradEqn, N, Q, U_x, U_y, U_z, mu, kappa, F)
         implicit none
         integer,          intent(in)  :: nEqn, nGradEqn
         integer         , intent(in)  :: N(3)
         real(kind=RP),    intent(in)  :: Q  (1:nEqn, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(in)  :: U_x(1:nGradEqn, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: U_y(1:nGradEqn, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: U_z(1:nGradEqn, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: mu  (0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(in)  :: kappa(0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(out) :: F   (1:nEqn, 0:N(1), 0:N(2), 0:N(3),1:NDIM)

         F(1,:,:,:,IX) = U_x(1,:,:,:)
         F(1,:,:,:,IY) = U_y(1,:,:,:)
         F(1,:,:,:,IZ) = U_z(1,:,:,:)

      end subroutine CHDivergenceFlux3D

      elemental subroutine AddQuarticDWPDerivative(c, mu)
         implicit none
         real(kind=RP), intent(in)    :: c
         real(kind=RP), intent(inout) :: mu
         
         mu = mu + 4.0_RP * c*(c-1.0_RP)*(c+1.0_RP) 

      end subroutine AddQuarticDWPDerivative

      elemental subroutine QuarticDWP(c, f)
         implicit none
         real(kind=RP), intent(in)  :: c
         real(kind=RP), intent(out) :: f

         f = POW2((c-1.0_RP)*(c+1.0_RP))

      end subroutine QuarticDWP

      pure subroutine PoiseuilleFlow(x, v)
         implicit none
         real(kind=RP), intent(in)  :: x(NDIM)
         real(kind=RP), intent(out) :: v(NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP), parameter  :: L = 40.0_RP

         v(1) = 0.0_RP
         v(2) = 0.0_RP
         v(3) = 4.0_RP * x(1) * ( L - x(1) )

      end subroutine PoiseuilleFlow

END Module Physics_CH
