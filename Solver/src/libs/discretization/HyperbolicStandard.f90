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