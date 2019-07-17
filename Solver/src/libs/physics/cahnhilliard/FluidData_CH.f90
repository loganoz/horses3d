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
      real(kind=RP)  :: eps      ! Coefficient in the dimensionless CH equation
      real(kind=RP)  :: invEps   ! (Inverse of the) Coefficient in the dimensionless CH equation
      real(kind=RP)  :: sigma    ! Interface energy
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
         ConstructMultiphase % eps          = 0.0_RP
         ConstructMultiphase % invEps       = 0.0_RP
         ConstructMultiphase % sigma        = 0.0_RP

      end function ConstructMultiphase
   
      subroutine SetMultiphase( multiphase_ )
         implicit none
         type(Multiphase_t), intent(in)  :: multiphase_

         multiphase % M       = multiphase_ % M
         multiphase % eps     = multiphase_ % eps
         multiphase % invEps  = 1.0_RP/multiphase % eps
         multiphase % sigma   = multiphase_ % sigma

      end subroutine SetMultiphase

end module FluidData_CH
