!
!//////////////////////////////////////////////////////
!
!   @File:    PhysicsStorage.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Apr 18 18:07:30 2018
!   @Last revision date: Thu Jul  5 12:34:56 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: feb27efbae31c25d40a6183082ebd1dcd742615e
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module PhysicsStorage
   use SMConstants, only: RP
#ifdef FLOW
   use FluidData, only: refValues, thermodynamics
#endif
#ifdef CAHNHILLIARD
   use FluidData, only: multiphase
#endif
#if defined(NAVIERSTOKES)
   use PhysicsStorage_NS
#elif defined(INCNS)
   use PhysicsStorage_iNS
#elif defined(MULTIPHASE)
   use PhysicsStorage_MU
#endif
#ifdef CAHNHILLIARD
   use PhysicsStorage_CH
#endif
   implicit none

#ifdef FLOW
   private refValues
#endif

#ifdef CAHNHILLIARD
   private multiphase
#endif
!
!  *****************************************************************************
!  These are the different modes supported by the ComputeTimeDerivative function
!        Defining this avoids to create specific procedures for each
!  *****************************************************************************
!
   enum, bind(C)
      enumerator :: CTD_IGNORE_MODE
#ifdef FLOW
      enumerator :: CTD_ONLY_NS
#ifdef CAHNHILLIARD
      enumerator :: CTD_NS_AND_CH
#endif
#endif
#ifdef CAHNHILLIARD
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
#elif defined(MULTIPHASE)
         call ConstructPhysicsStorage_MU( controlVariables, Lref, timeref, success )
#endif
!
!        Construct CHE physics
!        ---------------------
#ifdef CAHNHILLIARD
         call ConstructPhysicsStorage_CH( controlVariables, Lref, timeref, success )
#endif
!
!        ****************
!        Describe physics      
!        ****************
!
#if defined(NAVIERSTOKES)
         call DescribePhysicsStorage_NS()
#elif defined(INCNS)
         call DescribePhysicsStorage_iNS()
#elif defined(MULTIPHASE)
         call DescribePhysicsStorage_MU()
#endif

#ifdef CAHNHILLIARD
         call DescribePhysicsStorage_CH()
#endif

   
      end subroutine ConstructPhysicsStorage
end module PhysicsStorage
