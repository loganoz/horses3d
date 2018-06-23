!
!//////////////////////////////////////////////////////
!
!   @File:    HyperbolicStandard.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:16:32 2017
!   @Last revision date: Sat Jun 23 10:20:25 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: fce351220409e80ce5df1949249c2b870dd847aa
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
#if defined(NAVIERSTOKES) || defined(INCNS)
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
