!
!//////////////////////////////////////////////////////
!
!   @File:    FluidData.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Apr 18 18:07:28 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
module FluidData
#if defined(NAVIERSTOKES)
   use FluidData_NS
#endif

#if defined(CAHNHILLIARD)
   use FluidData_CH
#endif
   implicit none

end module FluidData
