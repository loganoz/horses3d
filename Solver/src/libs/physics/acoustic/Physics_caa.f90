!
!//////////////////////////////////////////////////////
!
!      CAA APE physics.
!      Modified from DSEM Code
!
!!     The variable mappings for the APE-4 Equations are
!!
!!              Q(1) = rho'
!!              Q(2) = u'
!!              Q(3) = v'
!!              Q(4) = w'
!!              Q(5) = p'
!
!////////////////////////////////////////////////////////////////////////
!    
!@mark -
!
#include "Includes.h"
!  **************
   module Physics_CAA
!  **************
!
      use SMConstants
      use PhysicsStorage_CAA
      use VariableConversion_CAA
      use FluidData_CAA
      use Utilities, only: outer_product
      implicit none

      private
      public  APEFlux, APEFlux1D
!
!     ========
      CONTAINS 
!     ========
!
!     
!
!//////////////////////////////////////////////////////////////////////////////
!
!           INVISCID FLUXES
!           ---------------   
!
!//////////////////////////////////////////////////////////////////////////////
!
!------------------------------------------------------------------------
!! Fluxes for APE I and IV. The density is not really solved, is computed
!! from pressure, this is why we have 0 flux.
!! for pressure: F = gamma*P_0*\vec{u'} + \vec{u_0}p'
!! for velocities: F = \vec{u_0} /dot /vec{u'} + p'/rho_0
!------------------------------------------------------------------------
!     
      pure subroutine APEFlux(Q, F, rho_, Qbase)
         implicit none
         real(kind=RP), intent(in)   :: Q(1:NCONS)
         real(kind=RP), intent(out)  :: F(1:NCONS , 1:NDIM)
         real(kind=RP), intent(in), optional :: rho_
         real(kind=RP), intent(in), optional :: Qbase(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
         real(kind=RP)              :: velocity_term

         velocity_term = dot_product( Qbase(ICAAU:ICAAW) , Q(ICAAU:ICAAW) ) + Q(ICAAP)/Q(ICAARHO)
!
!        X-Flux
!        ------         
         F(ICAARHO, IX ) = 0.0_RP
         F(ICAAU, IX ) = velocity_term
         F(ICAAV, IX ) = 0.0_RP
         F(ICAAW, IX ) = 0.0_RP
         F(ICAAP, IX ) = thermodynamics%gamma * Q(ICAAP) * Q(ICAAU) + Qbase(ICAAU) * Q(ICAAP)
!
!        Y-Flux
!        ------
         F(ICAARHO, IX ) = 0.0_RP
         F(ICAAU, IX ) = 0.0_RP
         F(ICAAV, IX ) = velocity_term
         F(ICAAW, IX ) = 0.0_RP
         F(ICAAP, IX ) = thermodynamics%gamma * Q(ICAAP) * Q(ICAAV) + Qbase(ICAAV) * Q(ICAAP)
!
!        Z-Flux
!        ------
         F(ICAARHO, IX ) = 0.0_RP
         F(ICAAU, IX ) = 0.0_RP
         F(ICAAV, IX ) = 0.0_RP
         F(ICAAW, IX ) = velocity_term
         F(ICAAP, IX ) = thermodynamics%gamma * Q(ICAAP) * Q(ICAAW) + Qbase(ICAAW) * Q(ICAAP)

      end subroutine APEFlux
!     
      pure subroutine APEFlux1D(Q, F, Qbase)
         implicit none
         real(kind=RP), intent(in)   :: Q(1:NCONS)
         real(kind=RP), intent(out)  :: F(1:NCONS)
         real(kind=RP), intent(in), optional :: Qbase(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------

!
!        only X-Flux
!        ------         
         F(ICAARHO, IX ) = 0.0_RP
         F(ICAAU, IX ) = dot_product( Qbase(ICAAU:ICAAW) , Q(ICAAU:ICAAW) ) + Q(ICAAP)/Q(ICAARHO)
         F(ICAAV, IX ) = 0.0_RP
         F(ICAAW, IX ) = 0.0_RP
         F(ICAAP, IX ) = thermodynamics%gamma * Q(ICAAP) * Q(ICAAU) + Qbase(ICAAU) * Q(ICAAP)

      end subroutine APEFlux1D
   END Module Physics_CAA
!
! /////////////////////////////////////////////////////////////////////
!
!----------------------------------------------------------------------
!! This routine returns the maximum eigenvalues for the Euler equations 
!! for the given solution value in each spatial direction. 
!! These are to be used to compute the local time step.
!----------------------------------------------------------------------
!
      SUBROUTINE ComputeEigenvaluesForStateCAA( Q, Qbase, eigen )
      
      USE SMConstants
      USE PhysicsStorage_CAA
      USE VariableConversion_CAA, ONLY:Pressure
      use FluidData_CAA,          only: Thermodynamics
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=Rp), DIMENSION(NCONS) :: Q
      REAL(KIND=Rp), DIMENSION(NCONS) :: Qbase
      REAL(KIND=Rp), DIMENSION(3)     :: eigen
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=Rp) :: u, v, w, p, a
!      
      associate ( gamma => thermodynamics % gamma ) 

      u = ABS( Q(ICAAU) )
      v = ABS( Q(ICAAV) )
      w = ABS( Q(ICAAW) )
      p = Pressure(Qbase)
      a = SQRT(gamma*p/Qbase(ICAARHO))
      
      eigen(1) = u + a
      eigen(2) = v + a
      eigen(3) = w + a

      end associate
      
      END SUBROUTINE ComputeEigenvaluesForState
