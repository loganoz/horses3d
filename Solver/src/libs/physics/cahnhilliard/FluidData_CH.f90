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
      real(kind=RP)  :: tCH           ! Chemical characteristic time
      real(kind=RP)  :: eps           ! Interface width
      real(kind=RP)  :: sigma         ! Interface tension
      real(kind=RP)  :: M0            ! Mobility
      real(kind=RP)  :: M0_star       ! Actual mobility
      real(kind=RP)  :: tCH_wDim      
      real(kind=RP)  :: eps_wDim    
      real(kind=RP)  :: sigma_wDim
      real(kind=RP)  :: M0_wDim
      real(kind=RP)  :: invEps        ! (Inverse of the) interface width
      contains
         procedure, nopass   :: SetStarMobility => Multiphase_SetStarMobility
   end type Multiphase_t

   type(Multiphase_t), protected    :: multiphase

   interface Multiphase_t
      module procedure ConstructMultiphase
   end interface Multiphase_t

   contains
      function ConstructMultiphase()
         implicit none
         type(Multiphase_t) :: ConstructMultiphase

         ConstructMultiphase % tCH        = 0.0_RP
         ConstructMultiphase % eps        = 0.0_RP
         ConstructMultiphase % sigma      = 0.0_RP
         ConstructMultiphase % M0         = 0.0_RP
         ConstructMultiphase % M0_star    = 0.0_RP
         ConstructMultiphase % tCH_wDim   = 0.0_RP
         ConstructMultiphase % eps_wDim   = 0.0_RP
         ConstructMultiphase % sigma_wDim = 0.0_RP
         ConstructMultiphase % M0_wDim    = 0.0_RP
         ConstructMultiphase % invEps     = 0.0_RP

      end function ConstructMultiphase
   
      subroutine SetMultiphase( multiphase_ )
         implicit none
         type(Multiphase_t), intent(in)  :: multiphase_

         multiphase % tCH        = multiphase_ % tCH
         multiphase % eps        = multiphase_ % eps
         multiphase % sigma      = multiphase_ % sigma
         multiphase % M0         = multiphase_ % M0
         multiphase % M0_star    = multiphase_ % M0
         multiphase % tCH_wDim   = multiphase_ % tCH_wDim
         multiphase % eps_wDim   = multiphase_ % eps_wDim
         multiphase % sigma_wDim = multiphase_ % sigma_wDim
         multiphase % M0_wDim    = multiphase_ % M0_wDim
         multiphase % invEps     = 1.0_RP / multiphase % eps

      end subroutine SetMultiphase
   
      subroutine Multiphase_SetStarMobility(M0)
         implicit none
         real(kind=RP),       intent(in) :: M0

         multiphase % M0_star = M0

      end subroutine Multiphase_SetStarMobility
end module FluidData_CH