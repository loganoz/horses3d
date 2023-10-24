#include "Includes.h"
module Physics_CH 
!
      USE SMConstants
      USE PhysicsStorage_CH
      use FluidData_CH, only: multiphase
      IMPLICIT NONE

      private
      public  CHDivergenceFlux, AddQuarticDWPDerivative, QuarticDWP
      public  Multiphase_AddChemFEDerivative
      public  PoiseuilleFlow
!
!     ---------
!     Constants
!     ---------
!
      INTEGER, PARAMETER   :: WALL_BC = 1, RADIATION_BC = 2
      INTEGER              :: boundaryCondition(4), bcType

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
      pure subroutine CHDivergenceFlux(nEqn, nGradEqn, Q, U_x, U_y, U_z, mu, beta, kappa, F)
         implicit none
         integer,       intent(in)  :: nEqn, nGradEqn
         real(kind=RP), intent(in)  :: Q   (nEqn)
         real(kind=RP), intent(in)  :: U_x (nGradEqn)
         real(kind=RP), intent(in)  :: U_y (nGradEqn)
         real(kind=RP), intent(in)  :: U_z (nGradEqn)
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: beta
         real(kind=RP), intent(in)  :: kappa
         real(kind=RP), intent(out) :: F(1:nEqn, 1:NDIM)

         F(1,IX) = U_x(1)
         F(1,IY) = U_y(1)
         F(1,IZ) = U_z(1)

      end subroutine CHDivergenceFlux

      elemental subroutine AddQuarticDWPDerivative(c, mu)
         implicit none
         real(kind=RP), intent(in)    :: c
         real(kind=RP), intent(inout) :: mu
         
         mu = mu + 4.0_RP * c*(c-1.0_RP)*(c+1.0_RP) 

      end subroutine AddQuarticDWPDerivative
   
      elemental subroutine Multiphase_AddChemFEDerivative(c, mu)
         implicit none
         real(kind=RP), intent(in)    :: c
         real(kind=RP), intent(inout) :: mu
         
         mu = mu + 48.0_RP * multiphase % sigma * multiphase % invEps * c*(c-1.0_RP)*(c-0.5_RP)

      end subroutine Multiphase_AddChemFEDerivative

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