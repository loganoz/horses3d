module ZoneClass
   use SMConstants
   use FaceClass
   use SharedBCModule
   use FTValueDictionaryClass
   use FTLinkedListClass
   
   
   private
   public Zone_t , ConstructZones, constructZoneModule

   integer, parameter      :: STR_LEN_ZONE = 128
   
   TYPE FTLinkedListPtr
      type(FTLinkedList), POINTER :: list
   END TYPE FTLinkedListPtr
   
   type Zone_t
      integer                     :: marker
      character(len=STR_LEN_ZONE) :: Name
      integer                     :: no_of_faces
      integer, allocatable        :: faces(:)
      contains
         procedure   :: Initialize => Zone_Initialize
   end type Zone_t
   
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
         use Headers
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
         character(len=STR_LEN_ZONE), allocatable :: zoneNames(:)
!
!        Get the number of markers from the Boundary Conditions dictionary
!        -----------------------------------------------------------------         
         no_of_markers = bcTypeDictionary % COUNT() 
         if ( no_of_markers .le. 0 ) return
!
!        Gather the zone names
!        ---------------------
         allocate ( zoneNames( 1:no_of_markers ) ) 
         zoneNames = bcTypeDictionary % allKeys()
!
!        Construct zones
!        ---------------
         allocate ( zones ( no_of_markers ) )
         do zoneID = 1 , no_of_markers
            call zones(zoneID) % Initialize ( zoneID , zoneNames(zoneID) )
         end do
!
!        Assign the faces
!        ----------------         
         call Zone_AssignFaces(faces,zones,no_of_markers,zoneNames)
         
         write(STD_OUT,'(/)')
         call Section_Header("Creating zones")
         write(STD_OUT,'(/)')
         
         do zoneID = 1, no_of_markers
            WRITE(STD_OUT,'(30X,A,A7,I0,A15,A)') "->", "  Zone ",zoneID, " for boundary: ",trim(zones(zoneID) % Name)
            write(STD_OUT,'(32X,A28,I0)') 'Number of faces: ', zones(zoneID) % no_of_faces
         end do
         
      end subroutine ConstructZones
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
         end do
!
!        Dump the linked list contents onto fixed-size arrays 
!        ----------------------------------------------------
         do zoneID = 1 , no_of_markers
            call zoneList(zoneID) % Load(realArray)
            allocate(zones(zoneID) % faces( size(realArray) ))
            zones(zoneID) % faces = floor(realArray)
            zones(zoneID) % no_of_faces = size(realArray)
            deallocate(realArray)
         end do
!
!        TODO: Destruct linked list
         
      end subroutine Zone_AssignFaces

end module ZoneClass
