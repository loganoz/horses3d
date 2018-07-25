!
!//////////////////////////////////////////////////////
!
!   @File:    EllipticDiscretizations.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:32:09 2017
!   @Last revision date: Sat Jun 23 10:20:23 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: fce351220409e80ce5df1949249c2b870dd847aa
!
!//////////////////////////////////////////////////////
!
module EllipticDiscretizations
   use EllipticDiscretizationClass
   use EllipticBR1
   use EllipticIP
   use EllipticBR2
   implicit none

   private
   public   EllipticFlux0D_f, EllipticFlux3D_f
   public   EllipticDiscretization_t, BassiRebay1_t, BassiRebay2_t, InteriorPenalty_t
   public   SIPG, IIPG, NIPG
#if defined(NAVIERSTOKES) || defined(INCNS)
   public   ViscousDiscretization
#endif
#if defined(CAHNHILLIARD)
   public   CHDiscretization
#endif


#if defined(NAVIERSTOKES) || defined(INCNS)
   class(EllipticDiscretization_t), allocatable  :: ViscousDiscretization
#endif
#if defined(CAHNHILLIARD)
   class(EllipticDiscretization_t), allocatable  :: CHDiscretization
#endif
      
end module EllipticDiscretizations
