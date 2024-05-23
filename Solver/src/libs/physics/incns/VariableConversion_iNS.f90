#include "Includes.h"
module VariableConversion_iNS
   use SMConstants
   use PhysicsStorage_iNS
   use FluidData_iNS
   implicit none

   private
   public   Pressure
   public   iNSGradientVariables
   public   GetiNSTwoFluidsViscosity, GetiNSOneFluidViscosity
   public   getVelocityGradients

   contains
!
! /////////////////////////////////////////////////////////////////////
!
!@mark -
!---------------------------------------------------------------------
!! Compute the pressure from the state variables
!---------------------------------------------------------------------
!
   PURE function Pressure(Q) RESULT(P)
   !
   !     ---------
   !     Arguments
   !     ---------
   !
         REAL(KIND=RP), DIMENSION(NCONS), INTENT(IN) :: Q
   !
   !     ---------------
   !     Local Variables
   !     ---------------
   !
         REAL(KIND=RP) :: P
         
         P = Q(INSP) 
   
         end function Pressure
   !
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! GradientValuesForQ takes the solution (Q) values and returns the
!! quantities of which the gradients will be taken.
!---------------------------------------------------------------------
!
      pure subroutine iNSGradientVariables(nEqn, nGrad, Q, U, rho_ )
         implicit none
         integer, intent(in)        :: nEqn, nGrad
         real(kind=RP), intent(in)  :: Q(nEqn)
         real(kind=RP), intent(out) :: U(nGrad)
         real(kind=RP), intent(in), optional :: rho_
!
!        ---------------
!        Local Variables
!        ---------------
!     
         real(kind=RP)  :: rho, invrho

         rho = Q(INSRHO)

         invRho = 1.0_RP / rho
         U(INSRHOU) = Q(INSRHOU) * invRho
         U(INSRHOV) = Q(INSRHOV) * invRho
         U(INSRHOW) = Q(INSRHOW) * invRho
         U(INSP)    = Q(INSP)
         U(INSRHO)  = -0.5_RP*(U(INSRHOU)*U(INSRHOU) + U(INSRHOV)*U(INSRHOV) + U(INSRHOW)*U(INSRHOW))

      end subroutine iNSGradientVariables

      pure subroutine GetiNSOneFluidViscosity(phi, mu)
!
!        ***********************************
!           Here phi is the density, such
!           that varies linearly from the
!           density of fluid 1 to that of
!           fluid 2
!        ***********************************
!
         implicit none
         real(kind=RP), intent(in)   :: phi
         real(kind=RP), intent(out)  :: mu
!
!        ---------------
!        Local variables
!        ---------------
!
         mu = dimensionless % mu(1)

      end subroutine GetiNSOneFluidViscosity

      pure subroutine GetiNSTwoFluidsViscosity(phi, mu)
!
!        ***********************************
!           Here phi is the density, such
!           that varies linearly from the
!           density of fluid 1 to that of
!           fluid 2
!        ***********************************
!
         implicit none
         real(kind=RP), intent(in)   :: phi
         real(kind=RP), intent(out)  :: mu
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)              :: mu2, mu1, rho1, rho2

         mu1 = dimensionless % mu(1)
         mu2 = dimensionless % mu(2)

         rho1 = dimensionless % rho(1)
         rho2 = dimensionless % rho(2)

         mu = mu1 * (phi - rho2)/(rho1-rho2) + mu2 * (phi-rho1)/(rho2-rho1)

      end subroutine GetiNSTwoFluidsViscosity



!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      pure subroutine getVelocityGradients(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
         implicit none
         !-arguments---------------------------------------------------
         real(kind=RP), intent(in)  :: Q(NCONS)
         real(kind=RP), intent(in)  :: Q_x(NGRAD), Q_y(NGRAD), Q_z(NGRAD)
         real(kind=RP), intent(out) :: U_x(NDIM), U_y(NDIM), U_z(NDIM)
         !-local-variables---------------------------------------------
         real(kind=RP) :: invRho, invRho2, uDivRho(NDIM)
         !-------------------------------------------------------------

         invRho  = 1._RP / Q(INSRHO)
         invRho2 = invRho * invRho
         
         uDivRho = [Q(INSRHOU) , Q(INSRHOV) , Q(INSRHOW) ] * invRho2
         
         u_x = invRho * Q_x(INSRHOU:INSRHOW) - uDivRho * Q_x(INSRHO)
         u_y = invRho * Q_y(INSRHOU:INSRHOW) - uDivRho * Q_y(INSRHO)
         u_z = invRho * Q_z(INSRHOU:INSRHOW) - uDivRho * Q_z(INSRHO)

      end subroutine getVelocityGradients
   
!
! /////////////////////////////////////////////////////////////////////
!
end module VariableConversion_iNS