#include "Includes.h"
module EllipticDiscretizations
   use EllipticDiscretizationClass
   use EllipticBR1
   use EllipticIP
   use EllipticBR2
   implicit none

   private
   public   EllipticFlux_f, GetViscosity_f
   public   EllipticDiscretization_t, BassiRebay1_t, BassiRebay2_t, InteriorPenalty_t
   public   SIPG, IIPG, NIPG
   public   ELLIPTIC_NS, ELLIPTIC_NSSA, ELLIPTIC_iNS, ELLIPTIC_CH, ELLIPTIC_MU
#ifdef FLOW
   public   ViscousDiscretization
#endif
#if defined(CAHNHILLIARD)
   public   CHDiscretization
#endif


#ifdef FLOW
   class(EllipticDiscretization_t), allocatable  :: ViscousDiscretization
#endif
#if defined(CAHNHILLIARD)
   class(EllipticDiscretization_t), allocatable  :: CHDiscretization
#endif
      
end module EllipticDiscretizations