!
!//////////////////////////////////////////////////////
!
!   @File:    HyperbolicStandard.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:16:32 2017
!   @Last revision date: Wed Jun 20 18:14:34 2018
!   @Last revision author: Juan Manzanero (j.manzanero1992@gmail.com)
!   @Last revision commit: 9c8ed8b6306ad0912cb55b510aa73d1610bb1cb5
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
