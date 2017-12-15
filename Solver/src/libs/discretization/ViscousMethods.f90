!
!//////////////////////////////////////////////////////
!
!   @File:    ViscousMethods.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:32:09 2017
!   @Last revision date: Fri Dec 15 10:18:28 2017
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 3eb47d0149dc29d08d7e725196514bb84fa71912
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
   public   ViscousMethod_t, BassiRebay1_t, InteriorPenalty_t, ViscousMethod
   public   BassiRebay2_t
!
   class(ViscousMethod_t), allocatable          :: ViscousMethod
!
!  ========
   contains
!  ========
!
end module ViscousMethods
