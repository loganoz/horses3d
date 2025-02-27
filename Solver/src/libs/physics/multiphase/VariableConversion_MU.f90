#include "Includes.h"
module VariableConversion_MU
   use SMConstants
   use PhysicsStorage_MU
   use FluidData_MU
   implicit none

   private
   public   mGradientVariables
   public   GetmTwoFluidsViscosity, GetmOneFluidViscosity
   public   getVelocityGradients

   contains
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! GradientValuesForQ takes the solution (Q) values and returns the
!! quantities of which the gradients will be taken.
!---------------------------------------------------------------------
!
      pure subroutine mGradientVariables( nEqn, nGrad, Q, U, rho_ )
!
!        --------------------------------------------------------------
!        Returns all gradient variables EXCEPT the chemical potential,
!        to be done manually afterwards.
!        --------------------------------------------------------------
!
         implicit none
         integer, intent(in)                 :: nEqn, nGrad
         real(kind=RP), intent(in)           :: Q(nEqn)
         real(kind=RP), intent(out)          :: U(nGrad)
         real(kind=RP), intent(in), optional :: rho_
!
!        ---------------
!        Local Variables
!        ---------------
!     
         real(kind=RP)  :: invSqrtRho

         invSqrtRho = 1.0_RP / sqrt(rho_)

         U = [0.0_RP, invSqrtRho, invSqrtRho, invSqrtRho, 1.0_RP] * Q

      end subroutine mGradientVariables

      pure subroutine GetmOneFluidViscosity(c, mu)
!
!        ***********************************
!           Here phi is the density, such
!           that varies linearly from the
!           density of fluid 1 to that of
!           fluid 2
!        ***********************************
!
         implicit none
         real(kind=RP), intent(in)   :: c
         real(kind=RP), intent(out)  :: mu
!
!        ---------------
!        Local variables
!        ---------------
!
         mu = dimensionless % mu(1)

      end subroutine GetmOneFluidViscosity

      pure subroutine GetmTwoFluidsViscosity(c, mu)
!
!        ***********************************
!           Here phi is the density, such
!           that varies linearly from the
!           density of fluid 1 to that of
!           fluid 2
!        ***********************************
!
         implicit none
         real(kind=RP), intent(in)   :: c
         real(kind=RP), intent(out)  :: mu
         real(kind=RP)  :: cHat

         cHat = min(max(c,0.0_RP),1.0_RP)

         mu = dimensionless % mu(1) * cHat + dimensionless % mu(2) * (1.0_RP - cHat)

      end subroutine GetmTwoFluidsViscosity
!
! /////////////////////////////////////////////////////////////////////
!
      pure subroutine getVelocityGradients(Q,Q_x,Q_y,Q_z,dimensionless_,U_x,U_y,U_z)
      implicit none
      !-arguments---------------------------------------------------
      real(kind=RP), intent(in)  :: Q(NCONS)
      real(kind=RP), intent(in)  :: Q_x(NGRAD), Q_y(NGRAD), Q_z(NGRAD)
      class(dimensionless_t),intent(in) :: dimensionless_ 
      real(kind=RP), intent(out) :: U_x(NDIM), U_y(NDIM), U_z(NDIM)
      !-local-variables---------------------------------------------
      real(kind=RP) :: rho, invRho, invRho2, uDivRho(NDIM)
      !-------------------------------------------------------------

       rho = dimensionless_ % rho(1)*Q(IMC) + dimensionless_ % rho(2)*(1-Q(IMC))
       invRho  = 1._RP / rho
       invRho2 = invRho * invRho
      
       uDivRho = [Q(IMSQRHOU) , Q(IMSQRHOV) , Q(IMSQRHOW) ] * invRho2
      
       u_x = invRho * Q_x(IMSQRHOU:IMSQRHOW) - uDivRho * Q_x(IMC)*(dimensionless_ % rho(1) - dimensionless_ % rho(2))
       u_y = invRho * Q_y(IMSQRHOU:IMSQRHOW) - uDivRho * Q_y(IMC)*(dimensionless_ % rho(1) - dimensionless_ % rho(2))
       u_z = invRho * Q_z(IMSQRHOU:IMSQRHOW) - uDivRho * Q_z(IMC)*(dimensionless_ % rho(1) - dimensionless_ % rho(2))

   end subroutine getVelocityGradients

!
! /////////////////////////////////////////////////////////////////////
!
end module VariableConversion_MU