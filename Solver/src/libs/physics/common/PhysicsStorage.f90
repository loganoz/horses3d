!
!//////////////////////////////////////////////////////
!
!   @File:    PhysicsStorage.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Apr 18 18:07:30 2018
!   @Last revision date: Sat Jun 23 10:20:32 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: fce351220409e80ce5df1949249c2b870dd847aa
!
!//////////////////////////////////////////////////////
!
module PhysicsStorage
   use SMConstants, only: RP
#if defined(NAVIERSTOKES)
   use PhysicsStorage_NS
#elif defined(INCNS)
   use PhysicsStorage_iNS
#endif
#if defined(CAHNHILLIARD)
   use PhysicsStorage_CH
#endif
   implicit none

#if (defined(NAVIERSTOKES) && !defined(CAHNHILLIARD))
   integer, parameter   :: NTOTALVARS = NCONS
   integer, parameter   :: NTOTALGRADS = NGRAD
#elif (!defined(NAVIERSTOKES) && defined(CAHNHILLIARD))
   integer, parameter   :: NTOTALVARS = NCOMP
   integer, parameter   :: NTOTALGRADS = NCOMP
#elif (defined(NAVIERSTOKES) && defined(CAHNHILLIARD))
   integer, parameter   :: NTOTALVARS = NCONS + NCOMP
   integer, parameter   :: NTOTALGRADS = NGRAD + NCOMP
#elif (defined(INCNS))
   integer, parameter   :: NTOTALVARS = NINC
   integer, parameter   :: NTOTALGRADS = NINC
#endif

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
#elif defined(INCNS)
         call ConstructPhysicsStorage_iNS( controlVariables, Lref, timeref, success )
#endif
!
!        Construct CHE physics
!        ---------------------
#if defined(CAHNHILLIARD)
         call ConstructPhysicsStorage_CH( controlVariables, Lref, timeref, success )
#endif
   
      end subroutine ConstructPhysicsStorage
end module PhysicsStorage
