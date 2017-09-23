module ZoneClass
   use SMConstants
   use FaceClass
   use SharedBCModule


   private
   public Zone_t , ConstructZones

   integer, parameter      :: STR_LEN_ZONE = 128

   type Zone_t
      integer                     :: marker
      character(len=STR_LEN_ZONE) :: Name
      integer                     :: no_of_faces
      integer, allocatable        :: faces(:)
      contains
         procedure   :: Initialize => Zone_Initialize
   end type Zone_t

   contains

      subroutine ConstructZones( faces , zones )
         implicit none
         class(Face), target                  :: faces(:)
         class(Zone_t), allocatable           :: zones(:)
!        --------------------------------------------------------
         integer                              :: zoneID
         class(FTMutableObjectArray), pointer :: bcObjects
         class(FTObject)            , pointer :: obj
         class(FTValue)             , pointer :: v
         integer                              :: no_of_markers

         bcObjects => bcTypeDictionary % allObjects()

         no_of_markers = bcObjects % COUNT()

         allocate ( zones ( 0 : no_of_markers ) )

         call zones(0) % Initialize( 0 , "Interior" )

         do zoneID = 1 , no_of_markers
            obj => bcObjects % objectAtIndex(zoneID)
            call castToValue(obj , v)
            call zones(zoneID) % Initialize ( zoneID , v % stringValue(requestedLength = STR_LEN_ZONE) )
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
