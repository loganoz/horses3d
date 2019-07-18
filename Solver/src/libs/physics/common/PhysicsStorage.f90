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
   use SMConstants, only: RP, STD_OUT
   use Headers
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

   real(kind=RP), protected     :: Lref
   real(kind=RP), protected     :: timeref

   character(len=*), parameter   :: REFERENCE_LENGTH_KEY = "reference length (m)" 

#if (!defined(FLOW)) && (defined(CAHNHILLIARD))
   integer, parameter :: NCONS = NCOMP
   integer, parameter :: NGRAD = NCOMP
#endif
   
   
   contains
      subroutine ConstructPhysicsStorage(controlVariables, success)
         USE FTValueDictionaryClass
         TYPE(FTValueDictionary)      :: controlVariables
         LOGICAL                      :: success
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: timeRef_NS, pRef
         

         if ( controlVariables % ContainsKey(REFERENCE_LENGTH_KEY) ) then
            Lref = controlVariables % DoublePrecisionValueForKey(REFERENCE_LENGTH_KEY)
         else
            Lref = 1.0_RP
         end if
!
!        Construct NSE physics
!        ---------------------
#if defined(NAVIERSTOKES)
         call ConstructPhysicsStorage_NS( controlVariables, Lref, timeRef_NS, success )
#elif defined(INCNS)
         call ConstructPhysicsStorage_iNS( controlVariables, Lref, timeRef_NS, success )
#elif defined(MULTIPHASE)
         call ConstructPhysicsStorage_MU( controlVariables, Lref, timeRef_NS, success )
#endif

!        Navier--Stokes equations set the reference time
!        -----------------------------------------------
#ifdef FLOW
         timeRef = timeRef_NS
         pRef    = refValues % p
#else
         timeRef = 1.0_RP
         pRef    = 1.0_RP
#endif
!
!        Construct CHE physics
!        ---------------------
#ifdef CAHNHILLIARD
         call ConstructPhysicsStorage_CH( controlVariables, Lref, timeRef, pRef, success )
#endif
!
!        ****************
!        Describe physics      
!        ****************
!
         call DescribePhysicsStorage_Common
#if defined(NAVIERSTOKES)
         call DescribePhysicsStorage_NS()
#elif defined(INCNS)
         call DescribePhysicsStorage_iNS()
#elif defined(MULTIPHASE)
         call DescribePhysicsStorage_MU()
#endif

#ifdef CAHNHILLIARD
         call DescribePhysicsStorage_CH(Lref)
#endif
   
      end subroutine ConstructPhysicsStorage

      subroutine DescribePhysicsStorage_Common()
         implicit none

         call Section_Header("Loading common physics")

         write(STD_OUT,'(/,/)')

         write(STD_OUT,'(30X,A,A40,ES10.3,A)') "->" , "Reference length (m): " , Lref
         write(STD_OUT,'(30X,A,A40,ES10.3,A)') "->" , "Reference time (s): "   , timeRef

      end subroutine DescribePhysicsStorage_Common

end module PhysicsStorage
