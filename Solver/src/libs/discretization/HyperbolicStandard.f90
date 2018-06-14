!
!//////////////////////////////////////////////////////
!
!   @File:    HyperbolicStandard.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:16:32 2017
!   @Last revision date: Wed Apr 18 20:19:03 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 0d746cd20d04ebda97f349d7f3b0b0fe00b5d7ca
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
#if defined(NAVIERSTOKES)
module HyperbolicStandard
   use SMConstants
   use RiemannSolvers_NS
   use HyperbolicDiscretizationClass
   implicit none

   private
   public   StandardDG_t

   type, extends(HyperbolicDiscretization_t)  :: StandardDG_t
   end type StandardDG_t

end module HyperbolicStandard
#endif
