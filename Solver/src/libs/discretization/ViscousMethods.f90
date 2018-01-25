!
!//////////////////////////////////////////////////////
!
!   @File:    ViscousMethods.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:32:09 2017
!   @Last revision date: Thu Jan 25 21:11:09 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 4ae0998f1881a7de77d8fb31fe8ac95dfed811ae
!
!//////////////////////////////////////////////////////
!
module ViscousMethods
   use ViscousMethodClass
   use ViscousBR1
   use ViscousIP
   use ViscousBR2
   implicit none
!
!
   private
   public   ViscousMethod_t, BassiRebay1_t, BassiRebay2_t, InteriorPenalty_t
   public   ViscousMethod,   BassiRebay1,   BassiRebay2,   InteriorPenalty
!
   class(ViscousMethod_t),  pointer :: ViscousMethod
   type(BassiRebay1_t),     target  :: BassiRebay1
   type(BassiRebay2_t),     target  :: BassiRebay2
   type(InteriorPenalty_t), target  :: InteriorPenalty
!
!  ========
   contains
!  ========
!
end module ViscousMethods
