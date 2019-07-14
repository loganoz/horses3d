!
!//////////////////////////////////////////////////////
!
!   @File:    HyperbolicStandard.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:16:32 2017
!   @Last revision date: Thu Sep 27 16:42:14 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 5ab4fc5764dead65069a92d809d881f964ea4900
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
#ifdef FLOW
module HyperbolicStandard
   use SMConstants
   use HyperbolicDiscretizationClass
   implicit none

   private
   public   StandardDG_t

   type, extends(HyperbolicDiscretization_t)  :: StandardDG_t
   end type StandardDG_t

end module HyperbolicStandard
#endif
