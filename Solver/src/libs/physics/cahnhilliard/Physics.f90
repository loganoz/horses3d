#include "Includes.h"
module Physics 
!
      USE SMConstants
      USE PhysicsStorage
      IMPLICIT NONE

      private
      public  ViscousFlux, QuarticDWPDerivative
!
!     ---------
!     Constants
!     ---------
!
      INTEGER, PARAMETER   :: WALL_BC = 1, RADIATION_BC = 2
      INTEGER              :: boundaryCondition(4), bcType

     interface ViscousFlux
      module procedure ViscousFlux0D, ViscousFlux2D, ViscousFlux3D
     end interface ViscousFlux

     interface QuarticDWPDerivative
      module procedure QuarticDWPDerivative_0D
      module procedure QuarticDWPDerivative_3D
     end interface QuarticDWPDerivative
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
      pure subroutine ViscousFlux0D(Q, U_x, U_y, U_z, mu, kappa, F)
         implicit none
         real(kind=RP), intent(in)  :: Q   (1:NCONS     )
         real(kind=RP), intent(in)  :: U_x (1:N_GRAD_EQN)
         real(kind=RP), intent(in)  :: U_y (1:N_GRAD_EQN)
         real(kind=RP), intent(in)  :: U_z (1:N_GRAD_EQN)
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: kappa
         real(kind=RP), intent(out) :: F(1:NCONS, 1:NDIM)

         F(1,IX) = U_x(1)
         F(1,IY) = U_y(1)
         F(1,IZ) = U_z(1)

      end subroutine ViscousFlux0D

      pure subroutine ViscousFlux2D( N, Q, U_x, U_y, U_z, mu, kappa, F)
         implicit none
         integer         , intent(in)  :: N(2)
         real(kind=RP),    intent(in)  :: Q  (1:NCONS, 0:N(1), 0:N(2))
         real(kind=RP),    intent(in)  :: U_x(1:N_GRAD_EQN, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: U_y(1:N_GRAD_EQN, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: U_z(1:N_GRAD_EQN, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: mu  (0:N(1), 0:N(2))
         real(kind=RP),    intent(in)  :: kappa(0:N(1), 0:N(2))
         real(kind=RP),    intent(out) :: F   (1:NCONS, 0:N(1), 0:N(2), 1:NDIM)

         F(1,:,:,IX) = U_x(1,:,:)
         F(1,:,:,IY) = U_y(1,:,:)
         F(1,:,:,IZ) = U_z(1,:,:)

      end subroutine ViscousFlux2D

      pure subroutine ViscousFlux3D( N, Q, U_x, U_y, U_z, mu, kappa, F)
         implicit none
         integer         , intent(in)  :: N(3)
         real(kind=RP),    intent(in)  :: Q  (1:NCONS, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(in)  :: U_x(1:N_GRAD_EQN, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: U_y(1:N_GRAD_EQN, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: U_z(1:N_GRAD_EQN, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: mu  (0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(in)  :: kappa(0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(out) :: F   (1:NCONS, 0:N(1), 0:N(2), 0:N(3),1:NDIM)

         F(1,:,:,:,IX) = U_x(1,:,:,:)
         F(1,:,:,:,IY) = U_y(1,:,:,:)
         F(1,:,:,:,IZ) = U_z(1,:,:,:)

      end subroutine ViscousFlux3D

      pure subroutine QuarticDWPDerivative_0D(c, f)
         implicit none
         real(kind=RP), intent(in)  :: c
         real(kind=RP), intent(out) :: f

         f = -c + POW3(c)

      end subroutine QuarticDWPDerivative_0D

      pure subroutine QuarticDWPDerivative_3D(N, c, mu)
         implicit none
         integer,       intent(in)    :: N(3)
         real(kind=RP), intent(in)    :: c (0:N(1),0:N(2),0:N(3))
         real(kind=RP), intent(inout) :: mu(0:N(1),0:N(2),0:N(3))

         mu = mu - c + POW3(c)

      end subroutine QuarticDWPDerivative_3D

END Module Physics
