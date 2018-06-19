!
!//////////////////////////////////////////////////////
!
!   @File:    BoundaryConditions.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Apr 18 18:07:28 2018
!   @Last revision date: Wed Jun 20 18:14:38 2018
!   @Last revision author: Juan Manzanero (j.manzanero1992@gmail.com)
!   @Last revision commit: 9c8ed8b6306ad0912cb55b510aa73d1610bb1cb5
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
