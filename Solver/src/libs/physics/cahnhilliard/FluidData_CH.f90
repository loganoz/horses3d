!
!//////////////////////////////////////////////////////
!
!   @File:    FluidData_CH.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Thu Apr 19 17:24:30 2018
!   @Last revision date: Thu Jul 26 17:26:20 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: ba557cd23630b1bd1f528599b9b33812f58d1f7b
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
      real(kind=RP)  :: eps      ! Coefficient in the dimensionless CH equation
      real(kind=RP)  :: invEps   ! (Inverse of the) Coefficient in the dimensionless CH equation
      real(kind=RP)  :: w        ! Dimensionless interface width
      real(kind=RP)  :: sigma    ! Interface energy
      real(kind=RP)  :: Pe       ! Peclet number (for CH)
      real(kind=RP)  :: Ca       ! Capilar number (for NS)
      real(kind=RP)  :: densityRatio 
      real(kind=RP)  :: viscRatio
      real(kind=RP)  :: barRho
      real(kind=RP)  :: tildeRho
      contains
         procedure :: SetDensityRatio   => Multiphase_SetDensityRatio
         procedure :: SetViscosityRatio => Multiphase_SetViscosityRatio
         procedure :: SetCapilarNumber  => Multiphase_SetCapilarNumber
   end type Multiphase_t

   type(Multiphase_t), protected    :: multiphase

   interface Multiphase_t
      module procedure ConstructMultiphase
   end interface Multiphase_t

   contains
      function ConstructMultiphase()
         implicit none
         type(Multiphase_t) :: ConstructMultiphase

         ConstructMultiphase % M            = 0.0_RP
         ConstructMultiphase % rhoS         = 0.0_RP
         ConstructMultiphase % kappa        = 0.0_RP
         ConstructMultiphase % c_alpha      = 0.0_RP
         ConstructMultiphase % c_beta       = 0.0_RP
         ConstructMultiphase % eps          = 0.0_RP
         ConstructMultiphase % invEps       = 0.0_RP
         ConstructMultiphase % w            = 0.0_RP
         ConstructMultiphase % sigma        = 0.0_RP
         ConstructMultiphase % Pe           = 0.0_RP
         ConstructMultiphase % Ca           = 0.0_RP
         ConstructMultiphase % densityRatio = 0.0_RP
         ConstructMultiphase % viscRatio    = 0.0_RP
         ConstructMultiphase % barRho       = 0.0_RP
         ConstructMultiphase % tildeRho     = 0.0_RP

      end function ConstructMultiphase
   
      subroutine SetMultiphase( multiphase_ )
         implicit none
         type(Multiphase_t), intent(in)  :: multiphase_

         multiphase % M       = multiphase_ % M
         multiphase % rhoS    = multiphase_ % rhoS
         multiphase % kappa   = multiphase_ % kappa
         multiphase % c_alpha = multiphase_ % c_alpha
         multiphase % c_beta  = multiphase_ % c_beta
         multiphase % eps     = multiphase_ % eps
         multiphase % invEps  = 1.0_RP/multiphase % eps
         multiphase % w       = multiphase_ % w
         multiphase % sigma   = multiphase_ % sigma
         multiphase % Pe      = multiphase_ % Pe
         multiphase % Ca      = multiphase_ % Ca
         multiphase % tildeRho  = multiphase_ % tildeRho
         multiphase % barRho    = multiphase_ % barRho

      end subroutine SetMultiphase

      subroutine Multiphase_SetDensityRatio(self, densityRatio)
         implicit none
         class(Multiphase_t)  :: self
         real(kind=RP), intent(in) :: densityRatio
   
         self % densityRatio = densityRatio

         self % tildeRho = 0.5_RP * (self % densityRatio - 1.0_RP)
         self % barRho   = 0.5_RP * (self % densityRatio + 1.0_RP)

      end subroutine Multiphase_SetDensityRatio

      subroutine Multiphase_SetViscosityRatio(self, viscRatio)
         implicit none
         class(Multiphase_t)  :: self
         real(kind=RP), intent(in) :: viscRatio

         self % viscRatio = viscRatio

      end subroutine Multiphase_SetViscosityRatio
   
      subroutine Multiphase_SetCapilarNumber(self, Ca)
         implicit none
         class(Multiphase_t)  :: self
         real(kind=RP), intent(in) :: Ca

         self % Ca = Ca

      end subroutine Multiphase_SetCapilarNumber
   
end module FluidData_CH
