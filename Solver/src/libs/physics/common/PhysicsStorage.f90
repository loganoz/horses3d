!
!//////////////////////////////////////////////////////
!
!   @File:    PhysicsStorage.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Apr 18 18:07:30 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
module PhysicsStorage
#if defined(NAVIERSTOKES)
   use PhysicsStorage_NS
#endif
#if defined(CAHNHILLIARD)
   use PhysicsStorage_CH
#endif
   implicit none

end module PhysicsStorage
