#include "Includes.h"
module BoundaryConditions
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
   implicit none

   private
   public   BCs, ConstructBoundaryConditions, DestructBoundaryConditions, SetBoundaryConditionsEqn
   public   NS_BC, C_BC, MU_BC

   enum, bind(C)
      enumerator :: INFLOW_BC = 1 , OUTFLOW_BC
      enumerator :: NOSLIPWALL_BC , FREESLIPWALL_BC
      enumerator :: PERIODIC_BC   , USERDEFINED_BC
   end enum

   character(len=BC_STRING_LENGTH), dimension(8)  :: implementedBCNames = [&
                "inflow              ",  &
                "outflow             ",  &
                "noslipwall          ",  &
                "freeslipwall        ",  &
                "periodic            ",  &
                "user-defined        ",  &
                "manufacturedsol     ",  &
                "msoutflowspecifyp   "]

   type BCSet_t
      class(GenericBC_t), allocatable :: bc
   end type BCSet_t

   type(BCSet_t), allocatable    :: BCs(:)
   integer                       :: no_of_BCs
   integer, protected            :: no_of_constructions = 0

   contains
      subroutine ConstructBoundaryConditions(no_of_zones, zoneNames)
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
!
!////////////////////////////////////////////////////////////////////////////
!
!        Procedures
!        ----------         
!
!////////////////////////////////////////////////////////////////////////////
!
      subroutine SetBoundaryConditionsEqn(which)
         implicit none
         integer,    intent(in) :: which
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: bcID

         do bcID = 1, no_of_BCs
            BCs(bcID) % bc % currentEqn = which
         end do

      end subroutine SetBoundaryConditionsEqn

      subroutine DestructBoundaryConditions()
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: bcID

         if ( no_of_constructions .gt. 1 ) then
            no_of_constructions = no_of_constructions - 1 
         else
            no_of_constructions = 0
            do bcID = 1, no_of_BCs
               call BCs(bcID) % bc % Destruct
            end do
      
            deallocate(BCs)
         end if

      end subroutine DestructBoundaryConditions
!
!////////////////////////////////////////////////////////////////////////////
!
!        Auxiliary subroutines
!        --------------------         
!
!////////////////////////////////////////////////////////////////////////////
!
      integer function GetZoneType(bname)
         implicit none
         character(len=*), intent(in)  :: bname
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: fid, io, bctype
         character(len=LINE_LENGTH) :: currentLine, loweredBname
         character(len=LINE_LENGTH) :: keyword, keyval
         logical                    :: inside
         type(FTValueDIctionary)    :: bcdict

         loweredbName = bname
         call toLower(loweredbName)
         call bcdict % initWithSize(16)

         open(newunit = fid, file = trim(controlFileName), status = "old", action = "read")
!
!        Navigate until the "#define boundary bname" sentinel is found
!        -------------------------------------------------------------
         inside = .false.
         do 
            read(fid, '(A)', iostat=io) currentLine

            IF(io .ne. 0 ) EXIT

            call PreprocessInputLine(currentLine)
            call toLower(currentLine)

            if ( index(trim(currentLine),"#define boundary") .ne. 0 ) then
               inside = CheckIfBoundaryNameIsContained(trim(currentLine), trim(loweredbname)) 
            end if
         
         
!
!           Get all keywords inside the zone
!           --------------------------------
            if ( inside ) then
               if ( trim(currentLine) .eq. "#end" ) exit

               keyword  = ADJUSTL(GetKeyword(currentLine))
               keyval   = ADJUSTL(GetValueAsString(currentLine))
               call ToLower(keyword)
      
               call bcdict % AddValueForKey(keyval, trim(keyword))

            end if

         end do

         if ( .not. bcdict % ContainsKey("type") ) then
            print*, "Missing boundary condition type for boundary ", trim(bname)
            errorMessage(STD_OUT)
            call exit(99)
         end if

         keyval = bcdict % StringValueForKey("type", LINE_LENGTH)

         call tolower(keyval)

         
         GetZoneType = -1
         do bctype = 1, size(implementedBCnames)
            if ( trim(keyval) .eq. trim(implementedBCnames(bctype)) ) then
               GetZoneType = bctype
            end if
         end do

         if ( GetZoneType .eq. -1 ) then
            print*, "Boundary type " ,trim(keyval), " not recognized."
            print*, "Options available are:"
            print*, "   * Inflow"
            print*, "   * Outflow"
            print*, "   * NoSlipWall"
            print*, "   * FreeSlipWall"
            print*, "   * Periodic"
            print*, "   * User-defined"
            errorMessage(STD_OUT)
            error stop
         end if

         call bcdict % Destruct
         close(fid)

      end function GetZoneType

end module BoundaryConditions