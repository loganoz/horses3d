!
!//////////////////////////////////////////////////////
!
!   @File:    PhysicsStorage.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Apr 18 18:07:30 2018
!   @Last revision date: Tue Jul  3 19:19:05 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 3db74c1b54d0c4fcf30b72bedefd8dbd2ef9b8ce
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
!
!  *****************************************************************************
!  These are the different modes supported by the ComputeTimeDerivative function
!        Defining this avoids to create specific procedures for each
!  *****************************************************************************
!
   enum, bind(C)
      enumerator :: CTD_IGNORE_MODE
#if defined(NAVIERSTOKES) || defined(INCNS)
      enumerator :: CTD_ONLY_NS
#if defined(CAHNHILLIARD)
      enumerator :: CTD_NS_AND_CH
#endif
#endif
#if defined(CAHNHILLIARD)
      enumerator :: CTD_ONLY_CH, CTD_ONLY_CH_LIN, CTD_ONLY_CH_NONLIN
#endif
   end enum

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
