!
!//////////////////////////////////////////////////////
!
!   @File:    BoundaryConditions.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Apr 18 18:07:28 2018
!   @Last revision date: Wed Apr 25 19:40:19 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 4749ed1216d5512d7b79f2485e9471f3161753ca
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

   enum, bind(C)
      enumerator :: NONE_BC
#if defined(NAVIERSTOKES)
      enumerator :: NS_BC
#endif
#if defined(CAHNHILLIARD)
      enumerator :: C_BC, MU_BC
#endif
   end enum

end module BoundaryConditionFunctions
