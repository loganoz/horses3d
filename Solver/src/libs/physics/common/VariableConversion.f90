!
!//////////////////////////////////////////////////////
!
!   @File:    VariableConversion.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Apr 18 18:07:30 2018
!   @Last revision date: Thu May 24 12:03:20 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: a9728294bcfa3ec9f4c553776074055792be41e2
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

   contains
      pure subroutine GetNSCHViscosity(phi, mu)
         use SMConstants, only: RP
         use FluidData
         implicit none
         real(kind=RP), intent(in)     :: phi
         real(kind=RP), intent(out)    :: mu
#if (defined(CAHNHILLIARD) && defined(NAVIERSTOKES))
         mu = 0.5_RP * dimensionless % mu * (1.0_RP - phi + multiphase % viscRatio * (1.0_RP + phi))
#else
         mu = 0.0_RP
#endif

      end subroutine GetNSCHViscosity

end module VariableConversion
