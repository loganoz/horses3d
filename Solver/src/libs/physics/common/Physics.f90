!
!//////////////////////////////////////////////////////
!
!   @File:    Physics.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Apr 18 18:07:29 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
module Physics
#if defined(NAVIERSTOKES)
   use Physics_NS
#endif
#if defined(CAHNHILLIARD)
   use Physics_CH
#endif
   implicit none

end module Physics
