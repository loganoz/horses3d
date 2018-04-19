!
!//////////////////////////////////////////////////////
!
!   @File:    PhysicsStorage.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Apr 18 18:07:30 2018
!   @Last revision date: Thu Apr 19 17:24:26 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: ca7f00098495d6fca03f13af3e8a139f88ed41e0
!
!//////////////////////////////////////////////////////
!
module PhysicsStorage
   use SMConstants, only: RP
#if defined(NAVIERSTOKES)
   use PhysicsStorage_NS
#endif
#if defined(CAHNHILLIARD)
   use PhysicsStorage_CH
#endif
   implicit none

   private RP

   real(kind=RP)     :: Lref
   real(kind=RP)     :: timeref
   
   contains
      subroutine ConstructPhysicsStorage(controlVariables, success)
         USE FTValueDictionaryClass
         TYPE(FTValueDictionary)      :: controlVariables
         LOGICAL                      :: success
!
!        Default values
!        --------------
         Lref = 1.0_RP
         timeref = 1.0_RP
!
!        Construct NSE physics
!        ---------------------
#if defined(NAVIERSTOKES)
         call ConstructPhysicsStorage_NS( controlVariables, Lref, timeref, success )
#endif
!
!        Construct CHE physics
!        ---------------------
#if defined(CAHNHILLIARD)
         call ConstructPhysicsStorage_CH( controlVariables, Lref, timeref, success )
#endif
   
      end subroutine ConstructPhysicsStorage
      
end module PhysicsStorage
