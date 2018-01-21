!
!//////////////////////////////////////////////////////
!
!   @File:    ViscousMethods.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:32:09 2017
!   @Last revision date: Sat Jan 20 18:43:02 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 1fb4c93edba89e57ab928b053c6b9b4598c17d82
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
