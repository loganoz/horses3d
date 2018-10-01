!
!//////////////////////////////////////////////////////
!
!   @File:    ZoneClass.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es), Andrés Rueda (am.rueda@upm.es)
!   @Created: 
!   @Last revision date: Mon Sep 24 19:27:04 2018
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 7ac2937102050656fd4a699d4a0b4a592b3431bf
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module ZoneClass
   use SMConstants
   use FaceClass
   use SharedBCModule
   use FTValueDictionaryClass
   use FTLinkedListClass
   use MeshTypes
   use BoundaryConditions, only: ConstructBoundaryConditions
   
   
   private
   public Zone_t , ConstructZones, ReassignZones, constructZoneModule, AllZoneNames

   integer, parameter      :: STR_LEN_ZONE = BC_STRING_LENGTH
   
   TYPE FTLinkedListPtr
      type(FTLinkedList), POINTER :: list
   END TYPE FTLinkedListPtr
   
   type Zone_t
      integer                     :: marker
      logical                     :: toBeDeleted = .false.
      character(len=STR_LEN_ZONE) :: Name
      integer                     :: no_of_faces
      integer, allocatable        :: faces(:)
      contains
         procedure   :: Initialize     => Zone_Initialize
         procedure   :: copy           => Zone_Assign
         generic     :: assignment(=)  => copy
   end type Zone_t
   
   contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ------------------------------------------
!     Creates the zones and assign faces to them 
!     according to values in the input file 
!     ------------------------------------------
      subroutine ConstructZones( faces , zones )
         use Headers
         use Utilities, only: toLower
         implicit none
         class(Face), target                  :: faces(:)
         class(Zone_t), allocatable           :: zones(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                  :: zoneID
         integer                                  :: no_of_markers
         character(len=STR_LEN_ZONE), pointer     :: zoneNames(:)
!
!        Get the number of markers from the Boundary Conditions dictionary
!        -----------------------------------------------------------------         
         no_of_markers = zoneNameDictionary % COUNT() 
         if ( no_of_markers .le. 0 ) return
!
!        Gather the zone names
!        ---------------------
         zoneNames => zoneNameDictionary % allKeys()
!
!        Construct zones
!        ---------------
         allocate ( zones ( no_of_markers ) )
         do zoneID = 1 , no_of_markers
            call toLower(zoneNames(zoneID))
            call zones(zoneID) % Initialize ( zoneID , zoneNames(zoneID) )
         end do
!
!        Assign the faces
!        ----------------         
         call Zone_AssignFaces(faces,zones,no_of_markers,zoneNames)
!
!        Construct Boundary conditions
!        ----------------------------- 
         call ConstructBoundaryConditions(no_of_markers, zoneNames)
         
         deallocate (zoneNames)
      end subroutine ConstructZones
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ------------------------------------------------------------------------------------------
!     Subroutine for reassigning the zone faces (needed when periodic boundaries are prescribed)
!     ------------------------------------------------------------------------------------------
      subroutine ReassignZones( faces , zones )
         use Headers
         implicit none
         class(Face), target                  :: faces(:)
         class(Zone_t)                        :: zones(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                  :: zoneID, fID
         integer                                  :: no_of_markers
         character(len=STR_LEN_ZONE), pointer     :: zoneNames(:)
!
!        Get the number of markers from the Boundary Conditions dictionary
!        -----------------------------------------------------------------         
         no_of_markers = zoneNameDictionary % COUNT() 
         if ( no_of_markers .le. 0 ) return
!
!        Gather the zone names
!        ---------------------
         zoneNames => zoneNameDictionary % allKeys()
!
!        Reset faces
!        -----------
         do fID = 1, size(faces)
            faces (fID) % zone = 0
         end do
!
!        Assign the faces
!        ----------------         
         call Zone_AssignFaces(faces,zones,no_of_markers,zoneNames)
         
         deallocate (zoneNames)
      end subroutine ReassignZones
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ------------------------------------------
!
!     ------------------------------------------
      subroutine Zone_Initialize ( self , marker , zoneName) 
         implicit none
         class(Zone_t)           :: self
         integer, intent(in)     :: marker
         character(len=*), intent(in)  :: zoneName

         self % marker = marker
         self % Name = trim(zoneName)
         self % no_of_faces = 0

      end subroutine Zone_Initialize
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------
!     Gather each zone faces
!     ----------------------
!
      subroutine Zone_AssignFaces(faces,zones,no_of_markers,zoneNames)
         use RealDataLinkedList
         implicit none
!        ------------------------------------------------
         integer                          :: no_of_markers
         type(Face)                       :: faces(:)
         type(Zone_t)                     :: zones(no_of_markers)
         character(len=STR_LEN_ZONE)      :: zoneNames(no_of_markers)
!        ------------------------------------------------
         integer                              :: fID, zoneID
         real(kind=RP), allocatable           :: realArray(:)
         type(RealDataLinkedList_t)   ,allocatable :: zoneList(:)
!
!        Initialize linked lists
!        -----------------------
         allocate(zoneList(no_of_markers))
!
!        Iterate over all faces and add the face ID to the zone linked list it belongs
!        -----------------------------------------------------------------------------
         do fID = 1, size(faces)
            if (faces(fID) % FaceType == HMESH_INTERIOR) cycle
            
            do zoneID = 1, no_of_markers
               if (trim(zoneNames(zoneID)) == trim(faces(fID) % boundaryName)) exit
            end do

            if (zoneID > no_of_markers) cycle
            
            call zoneList(zoneID) % Append( 1.0_RP * fID + 0.1_RP )
            faces (fID) % zone = zoneID
         end do
!
!        Dump the linked list contents onto fixed-size arrays 
!        ----------------------------------------------------
         do zoneID = 1 , no_of_markers
            call zoneList(zoneID) % Load(realArray)
            safedeallocate(zones(zoneID) % faces)
            allocate(zones(zoneID) % faces( size(realArray) ))
            zones(zoneID) % faces = floor(realArray)
            zones(zoneID) % no_of_faces = size(realArray)
            deallocate(realArray)
         end do
         
!        Finish up
!        ---------
         call zoneList % destruct
         deallocate (zoneList)
         
      end subroutine Zone_AssignFaces

      function AllZoneNames(no_of_zones, zones)
         implicit none
         integer,       intent(in)  :: no_of_zones
         class(Zone_t), intent(in)  :: zones(no_of_zones)
         character(len=LINE_LENGTH) :: AllZoneNames(no_of_zones)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer :: zID

         do zID = 1, no_of_zones
            AllZoneNames(zID) = trim(zones(zID) % Name)
         end do

      end function AllZoneNames
      
      elemental subroutine Zone_Assign (to, from)
         implicit none
         class(Zone_t), intent(inout)  :: to
         type(Zone_t) , intent(in)     :: from
         
         to % marker = from % marker
         to % toBeDeleted = from % toBeDeleted
         to % Name = from % Name
         to % no_of_faces = from % no_of_faces
         
         safedeallocate ( to % faces )
         allocate ( to % faces ( size(from % faces) ) )
         to % faces = from % faces
      end subroutine Zone_Assign
end module ZoneClass
