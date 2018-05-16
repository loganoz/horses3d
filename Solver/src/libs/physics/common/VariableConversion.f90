!
!//////////////////////////////////////////////////////
!
!   @File:    VariableConversion.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Apr 18 18:07:30 2018
!   @Last revision date: Fri Apr 20 17:25:10 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 056b1604b8f7d76486a7e001dc56e0b24c5e0edf
!
!//////////////////////////////////////////////////////
!
module VariableConversion
#if defined(NAVIERSTOKES)
   use VariableConversion_NS
#endif
#if defined(CAHNHILLIARD)
   use VariableConversion_CH
#endif
   implicit none

   abstract interface
      subroutine GetGradientValues0D_f(nEqn, nGradEqn, Q, U)
         use SMConstants, only: RP
         implicit none
         integer, intent(in)        :: nEqn, nGradEqn
         real(kind=RP), intent(in)  :: Q(nEqn)
         real(kind=RP), intent(out) :: U(nGradEqn)
      end subroutine GetGradientValues0D_f

      subroutine GetGradientValues3D_f(nEqn, nGradEqn, Nx, Ny, Nz, Q, U)
         use SMConstants, only: RP
         implicit none
         integer,    intent(in)  :: nEqn, nGradEqn, Nx, Ny, Nz
         real(kind=RP), intent(in)  :: Q(1:nEqn,  0:Nx, 0:Ny, 0:Nz)
         real(kind=RP), intent(out) :: U(1:nGradEqn, 0:Nx, 0:Ny, 0:Nz)
      end subroutine GetGradientValues3D_f
   end interface

end module VariableConversion
