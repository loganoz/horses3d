!
!//////////////////////////////////////////////////////
!
!   @File:    HyperbolicDiscretizations.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:16:31 2017
!   @Last revision date: Sat Jun 23 10:20:25 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: fce351220409e80ce5df1949249c2b870dd847aa
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
#if defined(NAVIERSTOKES) || defined(INCNS)
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
