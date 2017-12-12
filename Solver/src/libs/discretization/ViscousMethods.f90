!
!//////////////////////////////////////////////////////
!
!   @File:    ViscousMethods.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:32:09 2017
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
module ViscousMethods
   use ViscousMethodClass
   use ViscousBR1
   use ViscousIP
   implicit none
!
!
   private
   public   ViscousMethod_t , BassiRebay1_t , InteriorPenalty_t, ViscousMethod
!
   class(ViscousMethod_t), allocatable          :: ViscousMethod
!
!  ========
   contains
!  ========
!
end module ViscousMethods
