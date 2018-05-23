!
!//////////////////////////////////////////////////////
!
!   @File:    FluidData_CH.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Thu Apr 19 17:24:30 2018
!   @Last revision date: Wed May 23 12:57:26 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 7fde177b098184b58177a3a163cefdfebe7af55f
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module FluidData_CH
   use SMConstants
   implicit none

   private
   public   Multiphase_t, multiphase, SetMultiphase

   integer, parameter   :: STR_LEN_FLUIDDATA = 128
!
!  ----------------
!  Type definitions
!  ----------------
!
   type Multiphase_t
      real(kind=RP)  :: M        ! Mobility
      real(kind=RP)  :: rhoS     ! Double-well function height
      real(kind=RP)  :: kappa    ! Gradient energy coefficient
      real(kind=RP)  :: c_alpha  ! Alpha equilibrium concentration
      real(kind=RP)  :: c_beta   ! Beta equilibrium concentration
      real(kind=RP)  :: thetaw   ! Wall angle
      real(kind=RP)  :: eps      ! Coefficient in the dimensionless CH equation
      real(kind=RP)  :: w        ! Dimensionless interface width: 7.071 sqrt(kappa/rhoS)/Lref
      real(kind=RP)  :: sigma    ! Interface energy: 0.01508 sqrt(kappa*rhoS)/Lref
      real(kind=RP)  :: Pe       ! Peclet number (for CH)
      real(kind=RP)  :: Ca       ! Capilar number (for NS)
      real(kind=RP)  :: visc_ratio ! Viscosity ratio: mu(c=1)/mu(c=-1)
   end type Multiphase_t

   type(Multiphase_t), protected    :: multiphase

   contains
      subroutine SetMultiphase( multiphase_ )
         implicit none
         type(Multiphase_t), intent(in)  :: multiphase_

         multiphase % M       = multiphase_ % M
         multiphase % rhoS    = multiphase_ % rhoS
         multiphase % kappa   = multiphase_ % kappa
         multiphase % c_alpha = multiphase_ % c_alpha
         multiphase % c_beta  = multiphase_ % c_beta
         multiphase % thetaw  = multiphase_ % thetaw
         multiphase % eps     = multiphase_ % eps
         multiphase % w       = multiphase_ % w
         multiphase % sigma   = multiphase_ % sigma
         multiphase % Pe      = multiphase_ % Pe
         multiphase % Ca      = multiphase_ % Ca
         multiphase % visc_ratio = multiphase_ % visc_ratio

      end subroutine SetMultiphase

end module FluidData_CH
