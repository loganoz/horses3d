!
!//////////////////////////////////////////////////////
!
!   @File:    HyperbolicDiscretizations.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:16:31 2017
!   @Last revision date: Wed Apr 11 15:08:08 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 376f4a4922e7a9120e6c373514196f48dd0439c8
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
#if defined(NAVIERSTOKES)
module HyperbolicDiscretizations
   use SMConstants
   use RiemannSolvers
   use HyperbolicDiscretizationClass
   use HyperbolicStandard
   use HyperbolicSplitForm
   implicit none

   private
   public   HyperbolicDiscretization_t , StandardDG_t , SplitDG_t, HyperbolicDiscretization
   
   class(HyperbolicDiscretization_t), allocatable         :: HyperbolicDiscretization

end module HyperbolicDiscretizations
#endif
