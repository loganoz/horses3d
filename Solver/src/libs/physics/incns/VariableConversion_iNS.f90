!
!//////////////////////////////////////////////////////
!
!   @File:    VariableConversion_iNS.f90
!   @Author:  Juan Manzanero (j.manzanero1992@gmail.com)
!   @Created: Tue Jun 19 17:39:27 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
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
module VariableConversion_iNS
   use SMConstants
   use PhysicsStorage_iNS
   use FluidData_iNS
   implicit none

   private
   public   iNSGradientValuesForQ
   public   iNSGradientValuesForQ_0D, iNSGradientValuesForQ_3D
   public   GetiNSViscosity

   interface iNSGradientValuesForQ
       module procedure iNSGradientValuesForQ_0D , iNSGradientValuesForQ_3D
   end interface iNSGradientValuesForQ

   contains
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! GradientValuesForQ takes the solution (Q) values and returns the
!! quantities of which the gradients will be taken.
!---------------------------------------------------------------------
!
      pure subroutine iNSGradientValuesForQ_0D( nEqn, nGrad, Q, U )
         implicit none
         integer, intent(in)        :: nEqn, nGrad
         real(kind=RP), intent(in)  :: Q(nEqn)
         real(kind=RP), intent(out) :: U(nGrad)
!
!        ---------------
!        Local Variables
!        ---------------
!     
         U = Q

      end subroutine iNSGradientValuesForQ_0D

      pure subroutine iNSGradientValuesForQ_3D( nEqn, nGrad, Nx, Ny, Nz, Q, U )
         implicit none
         integer,       intent(in)  :: nEqn, nGrad, Nx, Ny, Nz
         real(kind=RP), intent(in)  :: Q(1:nEqn,  0:Nx, 0:Ny, 0:Nz)
         real(kind=RP), intent(out) :: U(1:nGrad, 0:Nx, 0:Ny, 0:Nz)

         U = Q

      end subroutine iNSGradientValuesForQ_3D

      pure subroutine GetiNSViscosity(phi, mu)
         implicit none
         real(kind=RP), intent(in)   :: phi
         real(kind=RP), intent(out)  :: mu

         mu = dimensionless % mu

      end subroutine GetiNSViscosity
!
! /////////////////////////////////////////////////////////////////////
!
end module VariableConversion_iNS
