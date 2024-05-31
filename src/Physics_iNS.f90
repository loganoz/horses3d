#include "Includes.h"
!  **************
   module Physics_iNS
!  **************
!
      use SMConstants
      use PhysicsStorage_iNS
      use VariableConversion_iNS
      use FluidData_iNS
      implicit none

      private
      public  iEulerFlux, iViscousFlux, iEulerXFlux
      public  getStressTensor
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
      pure subroutine iEulerFlux(Q, F, rho_)
         implicit none
         real(kind=RP), intent(in)   :: Q(1:NCONS)
         real(kind=RP), intent(out)  :: F(1:NCONS, 1:NDIM)
         real(kind=RP), intent(in), optional :: rho_
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: rho, invRho

         rho = (Q(INSRHO))

         invRho = 1.0_RP / rho
!
!        X-Flux
!        ------         
         F(INSRHO , IX) = Q(INSRHOU)
         F(INSRHOU, IX) = invRho*Q(INSRHOU)*Q(INSRHOU) + Q(INSP)
         F(INSRHOV, IX) = invRho*Q(INSRHOU)*Q(INSRHOV)
         F(INSRHOW, IX) = invRho*Q(INSRHOU)*Q(INSRHOW)
         F(INSP   , IX) = thermodynamics % rho0c02 * invRho * Q(INSRHOU)
!
!        Y-Flux
!        ------
         F(INSRHO , IY) = Q(INSRHOV)
         F(INSRHOU, IY) = invRho*Q(INSRHOV)*Q(INSRHOU)
         F(INSRHOV, IY) = invRho*Q(INSRHOV)*Q(INSRHOV) + Q(INSP)
         F(INSRHOW, IY) = invRho*Q(INSRHOV)*Q(INSRHOW)
         F(INSP   , IY) = thermodynamics % rho0c02 * invRho * Q(INSRHOV)
!
!        Z-Flux
!        ------
         F(INSRHO ,IZ) = Q(INSRHOW)
         F(INSRHOU,IZ) = invRho*Q(INSRHOW)*Q(INSRHOU)
         F(INSRHOV,IZ) = invRho*Q(INSRHOW)*Q(INSRHOV)
         F(INSRHOW,IZ) = invRho*Q(INSRHOW)*Q(INSRHOW) + Q(INSP)
         F(INSP   ,IZ) = thermodynamics % rho0c02 * invRho * Q(INSRHOW)
      
      end subroutine iEulerFlux

      pure subroutine iEulerXFlux(Q, F)
         implicit none
         real(kind=RP), intent(in)   :: Q(1:NCONS)
         real(kind=RP), intent(out)  :: F(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: rho, invRho

         rho = (Q(INSRHO))
      
         invRho = 1.0_RP / rho
!
!        X-Flux
!        ------         
         F(INSRHO ) = Q(INSRHOU)
         F(INSRHOU) = invRho*Q(INSRHOU)*Q(INSRHOU) + Q(INSP)
         F(INSRHOV) = invRho*Q(INSRHOU)*Q(INSRHOV)
         F(INSRHOW) = invRho*Q(INSRHOU)*Q(INSRHOW)
         F(INSP   ) = thermodynamics % rho0c02 * invRho * Q(INSRHOU)
      
      end subroutine iEulerXFlux
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
!>        VISCOUS FLUXES
!         --------------
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
      pure subroutine iViscousFlux(nEqn, nGradEqn, Q, U_x, U_y, U_z, mu, beta, kappa, F)
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

         F(INSRHO,IX)  = 0.0_RP
         F(INSRHOU,IX) = 2.0_RP * mu * U_x(INSRHOU)
         F(INSRHOV,IX) = mu * (U_x(INSRHOV) + U_y(INSRHOU))
         F(INSRHOW,IX) = mu * (U_x(INSRHOW) + U_z(INSRHOU))
         F(INSP,IX) = 0.0_RP

         F(INSRHO,IY)  = 0.0_RP
         F(INSRHOU,IY) = F(INSRHOV,IX)
         F(INSRHOV,IY) = 2.0_RP * mu * U_y(INSRHOV)
         F(INSRHOW,IY) = mu * (U_y(INSRHOW) + U_z(INSRHOV))
         F(INSP,IY) = 0.0_RP

         F(INSRHO,IZ)  = 0.0_RP
         F(INSRHOU,IZ) = F(INSRHOW,IX)
         F(INSRHOV,IZ) = F(INSRHOW,IY)
         F(INSRHOW,IZ) = 2.0_RP * mu * U_z(INSRHOW)
         F(INSP,IZ) = 0.0_RP

      end subroutine iViscousFlux

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      pure subroutine getStressTensor(Q,Q_x,Q_y,Q_z,tau)
         implicit none
         real(kind=RP), intent(in)      :: Q   (1:NCONS         )
         real(kind=RP), intent(in)      :: Q_x (1:NGRAD    )
         real(kind=RP), intent(in)      :: Q_y (1:NGRAD    )
         real(kind=RP), intent(in)      :: Q_z (1:NGRAD    )
         real(kind=RP), intent(out)     :: tau (1:NDIM, 1:NDIM   )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: divV
         real(kind=RP) :: U_x(NDIM), U_y(NDIM), U_z(NDIM), invRho,invRho2,uDivRho(NDIM)

         associate ( mu0 => dimensionless % mu )

            invRho  = 1._RP / Q(INSRHO)
            invRho2 = invRho * invRho
            
            uDivRho = [Q(INSRHOU) , Q(INSRHOV) , Q(INSRHOW) ] * invRho2

            U_x = invRho * Q_x(INSRHOU:INSRHOW) - uDivRho * Q_x(INSRHO)
            U_y = invRho * Q_y(INSRHOU:INSRHOW) - uDivRho * Q_y(INSRHO)
            U_z = invRho * Q_z(INSRHOU:INSRHOW) - uDivRho * Q_z(INSRHO)
            
            divV = U_x(IX) + U_y(IY) + U_z(IZ)
            
            ! tau(IX,IX) = mu0(IX) * (2.0_RP * U_x(IX) - 2.0_RP/3.0_RP * divV)
            ! tau(IY,IX) = mu0(IY) * ( U_x(IY) + U_y(IX) ) 
            ! tau(IZ,IX) = mu0(IZ) * ( U_x(IZ) + U_z(IX) ) 
            ! tau(IX,IY) = tau(IY,IX)
            ! tau(IY,IY) = mu0(IY) * (2.0_RP * U_y(IY) - 2.0_RP/3.0_RP * divV)
            ! tau(IZ,IY) = mu0(IZ) * ( U_y(IZ) + U_z(IY) ) 
            ! tau(IX,IZ) = tau(IZ,IX)
            ! tau(IY,IZ) = tau(IZ,IY)
            ! tau(IZ,IZ) = mu0(IZ) * (2.0_RP * U_z(IZ) - 2.0_RP/3.0_RP * divV)

            tau(IX,IX) = mu0(IX) * (2.0_RP * U_x(IX))
            tau(IY,IX) = mu0(IX) * ( U_x(IY) + U_y(IX) ) 
            tau(IZ,IX) = mu0(IX) * ( U_x(IZ) + U_z(IX) ) 
            tau(IX,IY) = tau(IY,IX)
            tau(IY,IY) = mu0(IX) * (2.0_RP * U_y(IY))
            tau(IZ,IY) = mu0(IX) * ( U_y(IZ) + U_z(IY) ) 
            tau(IX,IZ) = tau(IZ,IX)
            tau(IY,IZ) = tau(IZ,IY)
            tau(IZ,IZ) = mu0(IX) * (2.0_RP * U_z(IZ))
            
         end associate

      end subroutine getStressTensor

   END Module Physics_iNS
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
      
      USE SMConstants
      USE PhysicsStorage_iNS
      use FluidData_iNS,          only: Thermodynamics
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=Rp), DIMENSION(NCONS) :: Q
      REAL(KIND=Rp), DIMENSION(3)     :: eigen
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=Rp) :: u, v, w, p, a
!      
      u = ABS( Q(INSRHOU)/Q(INSRHO) )
      v = ABS( Q(INSRHOV)/Q(INSRHO) )
      w = ABS( Q(INSRHOW)/Q(INSRHO) )
      a = sqrt(u*u+v*v+w*w + 4.0_RP * thermodynamics % rho0c02/Q(INSRHO))
      
      eigen(1) = 0.5_RP * (u + a)
      eigen(2) = 0.5_RP * (v + a)
      eigen(3) = 0.5_RP * (w + a)
      
      END SUBROUTINE ComputeEigenvaluesForState