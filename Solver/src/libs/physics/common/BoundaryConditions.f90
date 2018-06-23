!
!//////////////////////////////////////////////////////
!
!   @File:    BoundaryConditions.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Apr 18 18:07:28 2018
!   @Last revision date: Sat Jun 23 10:20:30 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: fce351220409e80ce5df1949249c2b870dd847aa
!
!//////////////////////////////////////////////////////
!
module BoundaryConditionFunctions
#if defined(NAVIERSTOKES)
   use BoundaryConditionFunctions_NS
#elif defined(INCNS)
   use BoundaryConditionFunctions_iNS
#endif
#if defined(CAHNHILLIARD)
   use BoundaryConditionFunctions_CH
#endif
   implicit none

   enum, bind(C)
      enumerator :: NONE_BC
#if (defined(NAVIERSTOKES) || defined(INCNS))
      enumerator :: NS_BC
#endif
#if defined(CAHNHILLIARD)
      enumerator :: C_BC, MU_BC
#endif
   end enum

end module BoundaryConditionFunctions
