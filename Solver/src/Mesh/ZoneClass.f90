module ZoneClass
   use SMConstants
   use FaceClass
   use SharedBCModule
   
   
   private
   public Zone_t , ConstructZones, zoneNameDictionary, constructZoneModule

   integer, parameter      :: STR_LEN_ZONE = 128

   type Zone_t
      integer                     :: marker
      character(len=STR_LEN_ZONE) :: Name
      integer                     :: no_of_faces
      integer, allocatable        :: faces(:)
      contains
         procedure   :: Initialize => Zone_Initialize
   end type Zone_t
   
   TYPE(FTValueDictionary) :: zoneNameDictionary
   
   contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!   
!     -----------------------------------------------
!     Constructs the dictionary for storing the zones 
!     -----------------------------------------------
      subroutine constructZoneModule  
         implicit none
         call zoneNameDictionary % initWithSize(8)
      end subroutine constructZoneModule
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ------------------------------------------
!     Creates the zones and assign faces to them 
!     according to values in the input file 
!     ------------------------------------------
      subroutine ConstructZones( faces , zones )
         implicit none
         class(Face), target                  :: faces(:)
         class(Zone_t), allocatable           :: zones(:)
!        --------------------------------------------------------
         integer                              :: zoneID
         integer                              :: no_of_markers
         character(len=32), allocatable       :: zoneNames(:)
!        --------------------------------------------------------
         
         no_of_markers = zoneNameDictionary % COUNT()
         if (no_of_markers == 0) return               ! Only create zones if specified by the user
         
         allocate(zoneNames(no_of_markers))
         zoneNames = zoneNameDictionary % allKeys()
         
         allocate ( zones ( no_of_markers ) )
         
         !! Initialize zones

         do zoneID = 1 , no_of_markers
            call zones(zoneID) % Initialize ( zoneID , zoneNames(zoneID) )
            print*, "Zone ",zoneID , " created for boundary " , trim(zones(zoneID) % Name)
         end do
         
      end subroutine ConstructZones

      subroutine Zone_Initialize ( self , marker , zoneName) 
         implicit none
         class(Zone_t)           :: self
         integer, intent(in)     :: marker
         character(len=*), intent(in)  :: zoneName

         self % marker = marker
         self % Name = trim(zoneName)
         self % no_of_faces = 0

      end subroutine Zone_Initialize

end module ZoneClass
