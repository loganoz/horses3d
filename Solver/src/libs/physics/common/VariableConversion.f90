!
!//////////////////////////////////////////////////////
!
!   @File:    VariableConversion.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Apr 18 18:07:30 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
module VariableConversion
#if defined(NAVIERSTOKES)
   use VariableConversion_NS
#endif
#if defined(CAHNHILLIARD)
   use VariableConversion_CH
#endif
   implicit none

end module VariableConversion
