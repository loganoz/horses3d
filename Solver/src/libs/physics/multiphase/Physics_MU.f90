#include "Includes.h"
!  **************
   module Physics_MU
!  **************
!
      use SMConstants
      use PhysicsStorage_MU
      use VariableConversion_MU
      use FluidData_MU
      implicit none

      private
      public  mEulerFlux, mViscousFlux, mEulerXFlux

      real(kind=RP), parameter :: sigma_P = 0.0_RP

!
!     ========
      CONTAINS 
!     ========
!
!//////////////////////////////////////////////////////////////////////////////
!
!           INVISCID FLUXES
!           ---------------   
!
!//////////////////////////////////////////////////////////////////////////////
!
      pure subroutine mEulerFlux(Q, F, rho_, Qbase)
         implicit none
         real(kind=RP), intent(in)   :: Q(1:NCONS)
         real(kind=RP), intent(out)  :: F(1:NCONS, 1:NDIM)
         real(kind=RP), intent(in), optional :: rho_
         real(kind=RP), intent(in), optional :: Qbase(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: invSqrtRho

         invSqrtRho = 1.0_RP / sqrt(rho_)
!
!        X-Flux
!        ------         
         F(IMC     , IX) = Q(IMC)*invSqrtRho*Q(IMSQRHOU)
         F(IMSQRHOU, IX) = 0.5_RP*Q(IMSQRHOU)*Q(IMSQRHOU) + Q(IMP)
         F(IMSQRHOV, IX) = 0.5_RP*Q(IMSQRHOU)*Q(IMSQRHOV)
         F(IMSQRHOW, IX) = 0.5_RP*Q(IMSQRHOU)*Q(IMSQRHOW)
         F(IMP     , IX) = 0.0_RP
!
!        Y-Flux
!        ------
         F(IMC     , IY) = Q(IMC)*invSqrtRho*Q(IMSQRHOV)
         F(IMSQRHOU, IY) = 0.5_RP*Q(IMSQRHOV)*Q(IMSQRHOU)
         F(IMSQRHOV, IY) = 0.5_RP*Q(IMSQRHOV)*Q(IMSQRHOV) + Q(IMP)
         F(IMSQRHOW, IY) = 0.5_RP*Q(IMSQRHOV)*Q(IMSQRHOW)
         F(IMP     , IY) = 0.0_RP
!
!        Z-Flux
!        ------
         F(IMC     , IZ) = Q(IMC)*invSqrtRho*Q(IMSQRHOW)
         F(IMSQRHOU, IZ) = 0.5_RP*Q(IMSQRHOW)*Q(IMSQRHOU)
         F(IMSQRHOV, IZ) = 0.5_RP*Q(IMSQRHOW)*Q(IMSQRHOV) 
         F(IMSQRHOW, IZ) = 0.5_RP*Q(IMSQRHOW)*Q(IMSQRHOW) + Q(IMP)
         F(IMP     , IZ) = 0.0_RP
      
      end subroutine mEulerFlux

      pure subroutine mEulerXFlux(Q, F, rho)
         implicit none
         real(kind=RP), intent(in)   :: Q(1:NCONS)
         real(kind=RP), intent(out)  :: F(1:NCONS)
         real(kind=RP), intent(in)   :: rho
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: invSqrtRho

         invSqrtRho = 1.0_RP / sqrt(rho)
!
!        X-Flux
!        ------         
         F(IMC)      = Q(IMC)*invSqrtRho*Q(IMSQRHOU)
         F(IMSQRHOU) = 0.5_RP*Q(IMSQRHOU)*Q(IMSQRHOU) + Q(IMP)
         F(IMSQRHOV) = 0.5_RP*Q(IMSQRHOU)*Q(IMSQRHOV)
         F(IMSQRHOW) = 0.5_RP*Q(IMSQRHOU)*Q(IMSQRHOW)
         F(IMP)      = 0.0_RP
      
      end subroutine mEulerXFlux
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
!>        VISCOUS FLUXES
!         --------------
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
      pure subroutine mViscousFlux(nEqn, nGradEqn, Q, U_x, U_y, U_z, mu, beta, kappa, F)
         implicit none
         integer,       intent(in)  :: nEqn
         integer,       intent(in)  :: nGradEqn
         real(kind=RP), intent(in)  :: Q   (1:nEqn     )
         real(kind=RP), intent(in)  :: U_x (1:nGradEqn)
         real(kind=RP), intent(in)  :: U_y (1:nGradEqn)
         real(kind=RP), intent(in)  :: U_z (1:nGradEqn)
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: beta
         real(kind=RP), intent(in)  :: kappa
         real(kind=RP), intent(out) :: F(1:nEqn, 1:NDIM)

         F(IGMU,IX)  = beta*U_x(IGMU)
         F(IGU,IX) = 2.0_RP * mu * U_x(IGU)
         F(IGV,IX) = mu * (U_x(IGV) + U_y(IGU))
         F(IGW,IX) = mu * (U_x(IGW) + U_z(IGU))
         F(IGP,IX) = sigma_P*U_x(IGP)

         F(IGMU,IY)  = beta*U_y(IGMU)
         F(IGU,IY) = F(IGV,IX)
         F(IGV,IY) = 2.0_RP * mu * U_y(IGV)
         F(IGW,IY) = mu * (U_y(IGW) + U_z(IGV))
         F(IGP,IY) = sigma_P*U_y(IGP)

         F(IGMU,IZ)  = beta*U_z(IGMU)
         F(IGU,IZ) = F(IGW,IX)
         F(IGV,IZ) = F(IGW,IY)
         F(IGW,IZ) = 2.0_RP * mu * U_z(IGW)
         F(IGP,IZ) = sigma_P*U_z(IGP)

      end subroutine mViscousFlux

   END Module Physics_MU
!@mark -
!
! /////////////////////////////////////////////////////////////////////
!
!----------------------------------------------------------------------
!! This routine returns the maximum eigenvalues for the Euler equations 
!! for the given solution value in each spatial direction. 
!! These are to be used to compute the local time step.
!----------------------------------------------------------------------
!
      SUBROUTINE ComputeEigenvaluesForState( Q, eigen )
      
	  use SMConstants
      use PhysicsStorage_MU
      use VariableConversion_MU
      use FluidData_MU
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(NCONS) :: Q
      REAL(KIND=RP), DIMENSION(3)     :: eigen
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: u, v, w, p, a, invSqrtRho, sqrtRho
	  
      sqrtRho = sqrt(dimensionless % rho(2) + (dimensionless % rho(1)-dimensionless % rho(2))* (min(max(Q(1),0.0_RP),1.0_RP)))
	  invSqrtRho = 1.0_RP / sqrtRho
	   
      u = ABS( Q(IMSQRHOU)*invSqrtRho )
      v = ABS( Q(IMSQRHOV)*invSqrtRho )
      w = ABS( Q(IMSQRHOW)*invSqrtRho )
	  a = sqrt(thermodynamics % c02(2)+( thermodynamics % c02(1)-thermodynamics % c02(2))* (min(max(Q(1),0.0_RP),1.0_RP)))/ refValues % V
      
      eigen(1) = (u + a)
      eigen(2) = (v + a)
      eigen(3) = (w + a)
      
      END SUBROUTINE ComputeEigenvaluesForState
