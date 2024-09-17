#include "Includes.h"

submodule (BoundaryConditions) constructbcs
   use SMConstants
   use FTValueDictionaryClass,        only: FTValueDictionary
   use FileReaders,                   only: controlFileName
   use FileReadingUtilities,          only: GetKeyword, GetValueAsString, PreprocessInputLine
   use GenericBoundaryConditionClass, only: GenericBC_t, NS_BC, C_BC, MU_BC, CheckIfBoundaryNameIsContained
   use InflowBCClass,                 only: InflowBC_t
   use OutflowBCClass,                only: OutflowBC_t
   use NoSlipWallBCClass,             only: NoSlipWallBC_t
   use FreeSlipWallBCClass,           only: FreeSlipWallBC_t
   use PeriodicBCClass,               only: PeriodicBC_t
   use UserDefinedBCClass,            only: UserDefinedBC_t
   use Utilities, only: toLower, almostEqual
   use ZoneClass
   implicit none

   private
   
   contains
      
      module subroutine ConstructBoundaryConditions(no_of_zones, zoneNames)
         implicit none
         integer, intent(in)  :: no_of_zones
         character(len=*), intent(in)  :: zoneNames(no_of_zones)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: zID, zType

         no_of_constructions = no_of_constructions + 1
         if ( allocated(BCs) ) then
            print*, "*** WARNING!: Boundary conditions were previously allocated"
            return
         end if

         allocate(BCs(no_of_zones))
         no_of_BCs = no_of_zones

         do zID = 1, no_of_zones
!
!           Get zone type
!           -------------
            zType = GetZoneType(trim(zoneNames(zID)))
!
!           Allocate the boundary condition and construct
!           ---------------------------------------------
            select case(zType)
            case(INFLOW_BC)
               allocate(InflowBC_t       :: BCs(zID) % bc)
               select type(bc => BCs(zID) % bc)
               type is (InflowBC_t)
                  bc = InflowBC_t(trim(zoneNames(zID)))
               end select
            case(OUTFLOW_BC)
               allocate(OutflowBC_t      :: BCs(zID) % bc)
               select type(bc => BCs(zID) % bc)
               type is (OutflowBC_t)
                  bc = OutflowBC_t(trim(zoneNames(zID)))
               end select
            case(NOSLIPWALL_BC)
               allocate(NoSlipWallBC_t   :: BCs(zID) % bc)
               select type(bc => BCs(zID) % bc)
               type is (NoSlipWallBC_t)
                  bc = NoSlipWallBC_t(trim(zoneNames(zID)))
               end select
            case(FREESLIPWALL_BC)
               allocate(FreeSlipWallBC_t :: BCs(zID) % bc)
               select type(bc => BCs(zID) % bc)
               type is (FreeSlipWallBC_t)
                  bc = FreeSlipWallBC_t(trim(zoneNames(zID)))
               end select
            case(PERIODIC_BC)
               allocate(PeriodicBC_t     :: BCs(zID) % bc)
               select type(bc => BCs(zID) % bc)
               type is (PeriodicBC_t)
                  bc = PeriodicBC_t(trim(zoneNames(zID)))
               end select
            case(USERDEFINED_BC)
               allocate(UserDefinedBC_t  :: BCs(zID) % bc)
               select type(bc => BCs(zID) % bc)
               type is (UserDefinedBC_t)
                  bc = UserDefinedBC_t(trim(zoneNames(zID)))
               end select
            case default
               print*, "Unrecognized BC option"
               errorMessage(STD_OUT)
               error stop 99
            end select
             
         end do

      end subroutine ConstructBoundaryConditions

end submodule constructbcs