!
!//////////////////////////////////////////////////////
!
!   @File:    FluidData.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Apr 18 18:07:28 2018
!   @Last revision date: Thu Apr 19 17:24:25 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: ca7f00098495d6fca03f13af3e8a139f88ed41e0
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
