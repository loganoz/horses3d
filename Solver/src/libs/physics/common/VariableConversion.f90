!
!//////////////////////////////////////////////////////
!
!   @File:    VariableConversion.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Apr 18 18:07:30 2018
!   @Last revision date: Thu Jul  5 12:35:00 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: feb27efbae31c25d40a6183082ebd1dcd742615e
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module VariableConversion
#if defined(NAVIERSTOKES)
   use VariableConversion_NS
#elif defined(INCNS)
   use VariableConversion_iNS
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

#if (defined(CAHNHILLIARD) && defined(NAVIERSTOKES))
      pure subroutine GetNSCHViscosity(phi, mu)
         use SMConstants, only: RP
         use FluidData
         implicit none
         real(kind=RP), intent(in)     :: phi
         real(kind=RP), intent(out)    :: mu
!
!        ---------------
!        Local variables         
!        ---------------
!
         real(kind=RP)  :: cIn01, p

         cIn01 = 0.5_RP * (phi + 1.0_RP)
         p = POW3(cIn01) * (6.0_RP * POW2(cIn01) - 15.0_RP * cIn01 + 10.0_RP)

         mu = dimensionless % mu * ( (1.0_RP - p) + (p)*multiphase % viscRatio)

      end subroutine GetNSCHViscosity
#endif

#if (defined(INCNS) && defined(CAHNHILLIARD))
      pure subroutine GetiNSCHViscosity(phi, mu)
         use SMConstants, only: RP
         use FluidData
         implicit none
         real(kind=RP), intent(in)     :: phi
         real(kind=RP), intent(out)    :: mu
!
!        ---------------
!        Local variables         
!        ---------------
!
         real(kind=RP)  :: cIn01, p

         cIn01 = 0.5_RP * (phi + 1.0_RP)
         p = POW3(cIn01) * (6.0_RP * POW2(cIn01) - 15.0_RP * cIn01 + 10.0_RP)

         mu = dimensionless % mu(1) * (1.0_RP - p) + (p)*dimensionless % mu(2)
      
      end subroutine GetiNSCHViscosity
#endif

end module VariableConversion
