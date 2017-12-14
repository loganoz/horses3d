!
!//////////////////////////////////////////////////////
!
!   @File:    InviscidMethods.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:16:31 2017
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module InviscidMethods
   use SMConstants
   use RiemannSolvers
   use InviscidMethodClass
   use InviscidStandard
   use InviscidSplitForm
   implicit none

   private
   public   InviscidMethod_t , StandardDG_t , SplitDG_t, InviscidMethod
   
   class(InviscidMethod_t), allocatable         :: InviscidMethod
!
!  ========
   contains
!  ========
!
end module InviscidMethods
