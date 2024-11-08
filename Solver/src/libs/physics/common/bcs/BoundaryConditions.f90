#include "Includes.h"
module BoundaryConditions
   use SMConstants
   use FTValueDictionaryClass,        only: FTValueDictionary
   use FileReaders,                   only: controlFileName
   use SharedBCModule,                only: zoneNameDictionary
   use FileReadingUtilities,          only: GetKeyword, GetValueAsString, PreprocessInputLine
   ! use GenericBoundaryConditionClass, only: GenericBC_t, NS_BC, C_BC, MU_BC, CheckIfBoundaryNameIsContained
   use GenericBoundaryConditionClass, only: GenericBC_t, NS_BC, C_BC, MU_BC
   use InflowBCClass,                 only: InflowBC_t
   use OutflowBCClass,                only: OutflowBC_t
   use NoSlipWallBCClass,             only: NoSlipWallBC_t
   use FreeSlipWallBCClass,           only: FreeSlipWallBC_t
   use PeriodicBCClass,               only: PeriodicBC_t
   use UserDefinedBCClass,            only: UserDefinedBC_t
   use Utilities, only: toLower, almostEqual
   use ZoneClass, only: GetZoneType
   use MPI_Process_Info
   use HexMeshClass
#ifdef _HAS_MPI_
      use mpi
#endif
   implicit none

   private
   public   BCs, ConstructBoundaryConditions, DestructBoundaryConditions, SetBoundaryConditionsEqn, DescribeBoundaryConditions
   public   NS_BC, C_BC, MU_BC

   type BCSet_t
      class(GenericBC_t), allocatable :: bc
   end type BCSet_t

   type(BCSet_t), allocatable    :: BCs(:)
   integer                       :: no_of_BCs
   integer, protected            :: no_of_constructions = 0

   integer, parameter      :: STR_LEN_ZONE = BC_STRING_LENGTH

   contains
      subroutine ConstructBoundaryConditions()
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: zID, zType
         integer                                  :: no_of_zones
         character(len=STR_LEN_ZONE), pointer     :: zoneNames(:)
!
!        Get the number of markers from the Boundary Conditions dictionary
!        -----------------------------------------------------------------         
         no_of_zones = zoneNameDictionary % COUNT() 
         if ( no_of_zones .le. 0 ) return
!
!        Gather the zone names
!        ---------------------
         zoneNames => zoneNameDictionary % allKeys()
         
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

            call BCs(zID) % bc % CreateDeviceData()
             
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

      subroutine DescribeBoundaryConditions(mesh)
         USE Headers
         IMPLICIT NONE
         !-arguments------------------------------------------
         CLASS(HexMesh)      :: mesh
         !-local-variables------------------------------------
         integer           :: ierr
         integer           :: zoneID
         integer           :: no_of_bdry_faces
         integer           :: no_of_faces
         integer, allocatable :: facesPerZone(:)
         character(len=LINE_LENGTH) :: str
         !----------------------------------------------------
   
         allocate ( facesPerZone(size(mesh % zones)) )
   
   !     Gather information
   !     ------------------
   
         if (  MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
            do zoneID = 1, size(mesh % zones)
               call mpi_reduce ( mesh % zones(zoneID) % no_of_faces, facesPerZone(zoneID) , 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
            end do
   
            no_of_bdry_faces = sum(facesPerZone)
            no_of_faces      = (6*mesh % no_of_allElements + no_of_bdry_faces)/2
#endif
         else
            do zoneID = 1, size(mesh % zones)
               facesPerZone(zoneID) = mesh % zones(zoneID) % no_of_faces
            end do
   
            no_of_bdry_faces = sum(facesPerZone)
            no_of_faces = size ( mesh % faces )
         end if
   
   
   !     Describe the mesh
   !     -----------------
   
         if ( .not. MPI_Process % isRoot ) return
   
   !     Describe the zones
   !     ------------------
         write(STD_OUT,'(/)')
         call Section_Header("Creating zones")
         write(STD_OUT,'(/)')
   
         do zoneID = 1, size(mesh % zones)
            write(str,'(A,I0,A,A)') "Zone ", zoneID, " for boundary: ",trim(mesh % zones(zoneID) % Name)
            call SubSection_Header(trim(str))
            write(STD_OUT,'(30X,A,A28,I0)') "->", ' Number of faces: ', facesPerZone(zoneID)
            call BCs(zoneID) % bc % Describe
            write(STD_OUT,'(/)')
         end do
   
   !     Finish up
   !     ---------
         deallocate ( facesPerZone )

      end subroutine DescribeBoundaryConditions
end module BoundaryConditions
