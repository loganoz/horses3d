!
!//////////////////////////////////////////////////////
!
!   @File:    InviscidStandard.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:16:32 2017
!   @Last revision date: Tue Jan 16 11:59:32 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: cbae0faa7686246cad4b300efae466eda61473cd
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
#if defined(NAVIERSTOKES)
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
#endif
