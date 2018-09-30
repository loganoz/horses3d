!
!//////////////////////////////////////////////////////
!
!   @File:    HyperbolicDiscretizations.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:16:31 2017
!   @Last revision date: Thu Sep 27 16:42:13 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 5ab4fc5764dead65069a92d809d881f964ea4900
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
#if defined(NAVIERSTOKES) || defined(INCNS)
module HyperbolicDiscretizations
   use SMConstants
   use HyperbolicDiscretizationClass
   use HyperbolicStandard
   use HyperbolicSplitForm
   implicit none

   private
   public   HyperbolicDiscretization_t , StandardDG_t , SplitDG_t, HyperbolicDiscretization
   
   class(HyperbolicDiscretization_t), allocatable         :: HyperbolicDiscretization

end module HyperbolicDiscretizations
#endif
