!
!//////////////////////////////////////////////////////
!
!   @File:    BoundaryConditions.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Apr 18 18:07:28 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
module BoundaryConditionFunctions
#if defined(NAVIERSTOKES)
   use BoundaryConditionFunctions_NS
#endif
#if defined(CAHNHILLIARD)
   use BoundaryConditionFunctions_CH
#endif
   implicit none

end module BoundaryConditionFunctions
