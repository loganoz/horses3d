!
!//////////////////////////////////////////////////////
!
!   @File:    InviscidStandard.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:16:32 2017
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module InviscidStandard
   use SMConstants
   use RiemannSolvers
   use InviscidMethodClass
   implicit none

   private
   public   StandardDG_t

   type, extends(InviscidMethod_t)  :: StandardDG_t
   end type StandardDG_t

end module InviscidStandard
