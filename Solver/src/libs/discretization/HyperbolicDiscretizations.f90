#include "Includes.h"
#if defined(NAVIERSTOKES) || defined(INCNS)
#define HAS_SPLIT_FORMS
#endif

#ifdef FLOW
module HyperbolicDiscretizations
   use SMConstants
   use HyperbolicDiscretizationClass
   use HyperbolicStandard
#ifdef HAS_SPLIT_FORMS
   use HyperbolicSplitForm
   use HyperbolicSubcell
#endif
   implicit none

   private
   public HyperbolicDiscretization_t, StandardDG_t, HyperbolicDiscretization
#ifdef HAS_SPLIT_FORMS
   public SplitDG_t, SubcellDG_t
#endif
   
   class(HyperbolicDiscretization_t), allocatable :: HyperbolicDiscretization

end module HyperbolicDiscretizations
#endif
