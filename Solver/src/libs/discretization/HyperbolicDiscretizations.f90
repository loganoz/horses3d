!
!//////////////////////////////////////////////////////
!
!   @File:    HyperbolicDiscretizations.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:16:31 2017
!   @Last revision date: Wed Apr 18 20:19:01 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 0d746cd20d04ebda97f349d7f3b0b0fe00b5d7ca
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
#if defined(NAVIERSTOKES)
module HyperbolicDiscretizations
   use SMConstants
   use RiemannSolvers_NS
   use HyperbolicDiscretizationClass
   use HyperbolicStandard
   use HyperbolicSplitForm
   implicit none

   private
   public   HyperbolicDiscretization_t , StandardDG_t , SplitDG_t, HyperbolicDiscretization
   
   class(HyperbolicDiscretization_t), allocatable         :: HyperbolicDiscretization

end module HyperbolicDiscretizations
#endif
