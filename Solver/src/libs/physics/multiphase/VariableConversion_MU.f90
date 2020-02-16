!
!//////////////////////////////////////////////////////
!
!   @File:    VariableConversion_MU.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Jun 19 17:39:27 2018
!   @Last revision date: Mon Jul  2 14:17:29 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 7af1f42fb2bc9ea3a0103412145f2a925b4fac5e
!
!//////////////////////////////////////////////////////
!
!
!//////////////////////////////////////////////////////
!
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module VariableConversion_MU
   use SMConstants
   use PhysicsStorage_MU
   use FluidData_MU
   implicit none

   private
   public   mGradientVariables
   public   GetmTwoFluidsViscosity, GetmOneFluidViscosity

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
end module VariableConversion_MU
