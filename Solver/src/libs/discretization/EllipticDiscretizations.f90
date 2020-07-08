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
   public   ELLIPTIC_NS, ELLIPTIC_iNS, ELLIPTIC_CH, ELLIPTIC_MU
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
