!
!//////////////////////////////////////////////////////
!
!   @File:    Physics.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Apr 18 18:07:29 2018
!   @Last revision date: Sat Jun 23 10:20:31 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: fce351220409e80ce5df1949249c2b870dd847aa
!
!//////////////////////////////////////////////////////
!
module Physics
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
   use Physics_NS
#elif defined(SPALARTALMARAS)
	use Physics_NSSA
#elif defined(INCNS)
   use Physics_iNS
#elif defined(MULTIPHASE)
   use Physics_MU
#endif
#if defined(CAHNHILLIARD)
   use Physics_CH
#endif
   implicit none

end module Physics
