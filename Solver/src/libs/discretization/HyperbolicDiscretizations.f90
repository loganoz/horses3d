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
#define HAS_SPLIT_FORMS
#endif

#if defined(NAVIERSTOKES) || defined(INCNS) || defined(MULTIPHASE)
module HyperbolicDiscretizations
   use SMConstants
   use HyperbolicDiscretizationClass
   use HyperbolicStandard
#ifdef HAS_SPLIT_FORMS
   use HyperbolicSplitForm
#endif
   implicit none

   private
   public HyperbolicDiscretization_t, StandardDG_t, HyperbolicDiscretization
#ifdef HAS_SPLIT_FORMS
   public SplitDG_t
#endif
   
   class(HyperbolicDiscretization_t), allocatable         :: HyperbolicDiscretization

end module HyperbolicDiscretizations
#endif
