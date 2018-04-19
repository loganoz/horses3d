!
!//////////////////////////////////////////////////////
!
!   @File:    EllipticDiscretizations.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:32:09 2017
!   @Last revision date: Thu Jan 25 21:11:09 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 4ae0998f1881a7de77d8fb31fe8ac95dfed811ae
!
!//////////////////////////////////////////////////////
!
module EllipticDiscretizations
   use EllipticDiscretizationClass
   use EllipticBR1
   use EllipticIP
   use EllipticBR2
   implicit none
!
!
   private
   public   EllipticFlux0D_f, EllipticFlux3D_f
   public   EllipticDiscretization_t, BassiRebay1_t, BassiRebay2_t, InteriorPenalty_t
   public   EllipticDiscretization,   BassiRebay1,   BassiRebay2,   InteriorPenalty
   public   SIPG, IIPG, NIPG
!
   class(EllipticDiscretization_t),  pointer :: EllipticDiscretization
   type(BassiRebay1_t),     target  :: BassiRebay1
   type(BassiRebay2_t),     target  :: BassiRebay2
   type(InteriorPenalty_t), target  :: InteriorPenalty

      
end module EllipticDiscretizations
