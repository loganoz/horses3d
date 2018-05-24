!
!//////////////////////////////////////////////////////
!
!   @File:    VariableConversion_CH.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Thu Apr 19 17:24:31 2018
!   @Last revision date: Fri Apr 20 17:25:09 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 056b1604b8f7d76486a7e001dc56e0b24c5e0edf
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
      pure subroutine CHGradientValuesForQ_0D( nEqn, nGrad, Q, U )
         implicit none
         integer, intent(in)        :: nEqn, nGrad
         real(kind=RP), intent(in)  :: Q(nEqn)
         real(kind=RP), intent(out) :: U(nGrad)

         U = Q

      end subroutine CHGradientValuesForQ_0D

      pure subroutine CHGradientValuesForQ_3D( nEqn, nGrad, Nx, Ny, Nz, Q, U )
         implicit none
         integer,       intent(in)  :: nEqn, nGrad, Nx, Ny, Nz
         real(kind=RP), intent(in)  :: Q(1:nEqn,  0:Nx, 0:Ny, 0:Nz)
         real(kind=RP), intent(out) :: U(1:nGrad, 0:Nx, 0:Ny, 0:Nz)

         U = Q

      end subroutine CHGradientValuesForQ_3D

end module VariableConversion_CH
