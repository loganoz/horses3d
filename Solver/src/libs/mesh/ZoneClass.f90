module ZoneClass
   use SMConstants
   use FaceClass
   use SharedBCModule
   use FTLinkedListClass
   
   
   private
   public Zone_t , ConstructZones, zoneNameDictionary, constructZoneModule

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
         use Headers
         implicit none
         class(Face), target                  :: faces(:)
         class(Zone_t), allocatable           :: zones(:)
!        --------------------------------------------------------
         integer                              :: zoneID
         integer                              :: no_of_markers
         character(len=STR_LEN_ZONE), allocatable       :: zoneNames(:)
!        --------------------------------------------------------
         
         no_of_markers = zoneNameDictionary % COUNT()
         if (no_of_markers == 0) return               ! Only create zones if specified by the user
         
         allocate(zoneNames(no_of_markers))
         zoneNames = zoneNameDictionary % allKeys()
         
         allocate ( zones ( no_of_markers ) )
         
         !! Initialize zones

         do zoneID = 1 , no_of_markers
            call zones(zoneID) % Initialize ( zoneID , zoneNames(zoneID) )
         end do
         
         !! Assign faces to each zone
         call Zone_AssignFaces(faces,zones,no_of_markers,zoneNames)
         
         ! DEBUG
         write(STD_OUT,'(/)')
         call Section_Header("Creating zones")
         write(STD_OUT,'(/)')
         
         do zoneID = 1, no_of_markers
            WRITE(STD_OUT,'(30X,A,A7,I6,A15,A)') "->", "  Zone ",zoneID, " for boundary: ",trim(zones(zoneID) % Name)
            write(STD_OUT,'(32X,A28,I6)') 'Number of faces: ', zones(zoneID) % no_of_faces
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
!     -------------------------------------------------------
!     Subroutine for assigning faces to each zone, so that 
!     the drag and lift computations can be performed quickly 
!     -------------------------------------------------------
      subroutine Zone_AssignFaces(faces,zones,no_of_markers,zoneNames)
         implicit none
!        ------------------------------------------------
         integer                          :: no_of_markers
         type(Face)                       :: faces(:)
         type(Zone_t)                     :: zones(no_of_markers)
         character(len=STR_LEN_ZONE)      :: zoneNames(no_of_markers)
!        ------------------------------------------------
         integer                              :: fID, zoneID
         type(FTLinkedListPtr)   ,allocatable :: zoneList(:)
         CLASS(FTValue)             , POINTER :: v
         CLASS(FTObject)            , POINTER :: objectPtr
         type(FTMutableObjectArray) , POINTER :: array
!        ------------------------------------------------
!
!        --------
!        Initialize linked lists
!        --------
!
         allocate(zoneList(no_of_markers))
         
         do zoneID = 1, no_of_markers
            allocate (zoneList(zoneID) % list)
            call zoneList(zoneID) % list % init()
         end do
!
!        --------
!        We first iterate over all faces and add the face ID to the corresponding linked list
!        --------
!
         do fID = 1, size(faces)
            if (faces(fID) % FaceType == HMESH_INTERIOR) cycle
            
            do zoneID = 1, no_of_markers
               if (trim(zoneNames(zoneID)) == trim(faces(fID) % boundaryName)) exit
            end do
            if (zoneID > no_of_markers) cycle
            
            allocate (v)
            call v % initWithValue(fID)
            objectPtr => v
            call zoneList(zoneID) % list % add (objectPtr)
            CALL release(v)
            
         end do
!
!        --------
!        Now we create the arrays with the information in the linked lists
!        --------
!
         do zoneID = 1 , no_of_markers
            array => zoneList(zoneID) % list % allObjects()
            
            zones(zoneID) % no_of_faces = array % COUNT()
            allocate (zones(zoneID) % faces (zones(zoneID) % no_of_faces))
            
            do fID = 1, zones(zoneID) % no_of_faces
               
               v => valueFromObject(array % objectAtIndex(fID))
               zones(zoneID) % faces(fID) = v % integerValue()
               
            end do
            
         end do
         
!
!        --------
!        Clean up
!        --------
!
         call release(array)
         nullify(objectPtr)
         nullify(v)
         
         do zoneID = 1, no_of_markers
            call release (zoneList(zoneID) % list)
         end do
         
      end subroutine Zone_AssignFaces

end module ZoneClass
