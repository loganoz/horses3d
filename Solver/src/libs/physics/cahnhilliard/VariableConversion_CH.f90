!
!//////////////////////////////////////////////////////
!
!   @File:    VariableConversion_CH.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Thu Apr 19 17:24:31 2018
!   @Last revision date: Thu May 24 12:03:20 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: a9728294bcfa3ec9f4c553776074055792be41e2
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
module VariableConversion_CH
   use SMConstants
   use PhysicsStorage_CH
   implicit none

   private
   public CHGradientValuesForQ, CHGradientValuesForQ_0D, CHGradientValuesForQ_3D
   public GetCHViscosity

   interface CHGradientValuesForQ
       module procedure CHGradientValuesForQ_0D , CHGradientValuesForQ_3D
   end interface CHGradientValuesForQ

   contains
!
!---------------------------------------------------------------------
!! GradientValuesForQ takes the solution (Q) values and returns the
!! quantities of which the gradients will be taken.
!---------------------------------------------------------------------
!
      pure subroutine CHGradientValuesForQ_0D( nEqn, nGrad, Q, U, rho_ )
         implicit none
         integer, intent(in)        :: nEqn, nGrad
         real(kind=RP), intent(in)  :: Q(nEqn)
         real(kind=RP), intent(out) :: U(nGrad)
         real(kind=RP), intent(in), optional :: rho_

         U = Q

      end subroutine CHGradientValuesForQ_0D

      pure subroutine CHGradientValuesForQ_3D( nEqn, nGrad, Nx, Ny, Nz, Q, U, rho_ )
         implicit none
         integer,       intent(in)  :: nEqn, nGrad, Nx, Ny, Nz
         real(kind=RP), intent(in)  :: Q(1:nEqn,  0:Nx, 0:Ny, 0:Nz)
         real(kind=RP), intent(out) :: U(1:nGrad, 0:Nx, 0:Ny, 0:Nz)
         real(kind=RP), intent(in), optional :: rho_(0:Nx, 0:Ny, 0:Nz)

         U = Q

      end subroutine CHGradientValuesForQ_3D

      pure subroutine GetCHViscosity(phi, mu)
         implicit none
         real(kind=RP), intent(in)  :: phi
         real(kind=RP), intent(out) :: mu

         mu = 1.0_RP

      end subroutine GetCHViscosity

end module VariableConversion_CH
