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
module PhysicsStorage
   use SMConstants, only: RP
#if (defined(NAVIERSTOKES) || defined(INCNS))
   use FluidData, only: refValues, thermodynamics
#endif
#if defined(CAHNHILLIARD)
   use FluidData, only: multiphase
#endif
#if defined(NAVIERSTOKES)
   use PhysicsStorage_NS
#elif defined(INCNS)
   use PhysicsStorage_iNS
#endif
#if defined(CAHNHILLIARD)
   use PhysicsStorage_CH
#endif
   implicit none

#if defined(NAVIERSTOKES) || defined(INCNS)
   private refValues
#endif

#if defined(CAHNHILLIARD)
   private multiphase
#endif

#if (defined(NAVIERSTOKES) && !defined(CAHNHILLIARD))
   integer, parameter   :: NTOTALVARS = NCONS
   integer, parameter   :: NTOTALGRADS = NGRAD
#elif (defined(NAVIERSTOKES) && defined(CAHNHILLIARD))
   integer, parameter   :: NTOTALVARS = NCONS + NCOMP
   integer, parameter   :: NTOTALGRADS = NGRAD + NCOMP
#elif (defined(INCNS) && !defined(CAHNHILLIARD))
   integer, parameter   :: NTOTALVARS = NCONS
   integer, parameter   :: NTOTALGRADS = NCONS
#elif (defined(INCNS) && defined(CAHNHILLIARD))
   integer, parameter   :: NTOTALVARS = NCONS + NCOMP
   integer, parameter   :: NTOTALGRADS = NCONS + NCOMP
#elif defined(CAHNHILLIARD)
   integer, parameter   :: NTOTALVARS = NCOMP
   integer, parameter   :: NTOTALGRADS = NCOMP
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

#if (defined(NAVIERSTOKES) && defined(CAHNHILLIARD))
   character(len=*), parameter   :: DENSITY_RATIO_KEY   = "density ratio (rho2/rho1)"
   character(len=*), parameter   :: VISCOSITY_RATIO_KEY = "viscosity ratio (mu2/mu1)"
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
!
!        Define density and viscosity ratios
!        -----------------------------------
#if (defined(NAVIERSTOKES) && defined(CAHNHILLIARD)) 
         if ( controlVariables % ContainsKey(DENSITY_RATIO_KEY) ) then
            call multiphase % SetDensityRatio(controlVariables % DoublePrecisionValueForKey(DENSITY_RATIO_KEY))
         else
            call multiphase % SetDensityRatio(1.0_RP)
         end if

         if ( controlVariables % ContainsKey(VISCOSITY_RATIO_KEY) ) then
            call multiphase % SetViscosityRatio(controlVariables % DoublePrecisionValueForKey(VISCOSITY_RATIO_KEY))
         else
            call multiphase % SetViscosityRatio(1.0_RP)
         end if
#endif

#if (defined(INCNS) && defined(CAHNHILLIARD)) 
         call multiphase % SetDensityRatio(thermodynamics % rho(2) / thermodynamics % rho(1))
         call multiphase % SetViscosityRatio(thermodynamics % mu(2)  / thermodynamics % mu(1))
#endif


#if (defined(NAVIERSTOKES) || defined(INCNS)) && defined(CAHNHILLIARD)
!
!        Compute the capilar number
!        --------------------------
         call multiphase % SetCapilarNumber( (2.0_RP * sqrt(2.0_RP) / 3.0_RP) * refValues % V * refValues % mu / multiphase % sigma)
         
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
#endif

#if defined(CAHNHILLIARD)
         call DescribePhysicsStorage_CH()
#endif

   
      end subroutine ConstructPhysicsStorage
end module PhysicsStorage
