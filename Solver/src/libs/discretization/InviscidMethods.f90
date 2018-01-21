!
!//////////////////////////////////////////////////////
!
!   @File:    InviscidMethods.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:16:31 2017
!   @Last revision date: Tue Jan 16 11:59:31 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: cbae0faa7686246cad4b300efae466eda61473cd
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
#if defined(NAVIERSTOKES)
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
#endif
