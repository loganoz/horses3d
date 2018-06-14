!
!//////////////////////////////////////////////////////
!
!   @File:    EllipticDiscretizations.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:32:09 2017
!   @Last revision date: Wed May 23 12:57:24 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 7fde177b098184b58177a3a163cefdfebe7af55f
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
#if defined(NAVIERSTOKES)
   public   ViscousDiscretization
#endif
#if defined(CAHNHILLIARD)
   public   CHDiscretization
#endif


#if defined(NAVIERSTOKES)
   class(EllipticDiscretization_t), allocatable  :: ViscousDiscretization
#endif
#if defined(CAHNHILLIARD)
   class(EllipticDiscretization_t), allocatable  :: CHDiscretization
#endif
      
end module EllipticDiscretizations
