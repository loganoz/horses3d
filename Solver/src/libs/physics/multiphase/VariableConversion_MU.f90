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
   public   mGradientValuesForQ
   public   mGradientValuesForQ_0D, mGradientValuesForQ_3D
   public   GetmTwoFluidsViscosity, GetmOneFluidViscosity

   interface mGradientValuesForQ
       module procedure mGradientValuesForQ_0D , mGradientValuesForQ_3D
   end interface mGradientValuesForQ

   contains
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! GradientValuesForQ takes the solution (Q) values and returns the
!! quantities of which the gradients will be taken.
!---------------------------------------------------------------------
!
      pure subroutine mGradientValuesForQ_0D( nEqn, nGrad, Q, U, rho_ )
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

         U(IGU) = Q(IMSQRHOU) * invSqrtRho
         U(IGV) = Q(IMSQRHOV) * invSqrtRho
         U(IGW) = Q(IMSQRHOW) * invSqrtRho
         U(IGP) = Q(IMP)

      end subroutine mGradientValuesForQ_0D

      pure subroutine mGradientValuesForQ_3D( nEqn, nGrad, Nx, Ny, Nz, Q, U, rho_ )
         implicit none
         integer,       intent(in)  :: nEqn, nGrad, Nx, Ny, Nz
         real(kind=RP), intent(in)  :: Q(1:nEqn,  0:Nx, 0:Ny, 0:Nz)
         real(kind=RP), intent(out) :: U(1:nGrad, 0:Nx, 0:Ny, 0:Nz)
         real(kind=RP), intent(in), optional :: rho_(0:Nx, 0:Ny, 0:Nz)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, k
         real(kind=RP) :: invSqrtRho

         do k = 0, Nz ; do j = 0, Ny ; do i = 0, Nx
            invSqrtRho = 1.0_RP / sqrt(rho_(i,j,k))
!
!           I made this an entire line just in case the compiler vectorizes it ?
!           ------------------------------------------------------------------
            U(IGU:IGP,i,j,k) = [invSqrtRho, invSqrtRho, invSqrtRho, 1.0_RP] * Q(IMSQRHOU:IMP,i,j,k)
            
         end do       ; end do       ; end do
         
      end subroutine mGradientValuesForQ_3D

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

         mu = dimensionless % mu(1) * c + dimensionless % mu(2) * (1.0_RP - c)

      end subroutine GetmTwoFluidsViscosity

!
! /////////////////////////////////////////////////////////////////////
!
end module VariableConversion_MU
