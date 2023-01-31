!
!////////////////////////////////////////////////////////////////////////
!
!      FTDictionary.f90
!      Created: January 28, 2013 2:00 PM 
!      By: David Kopriva 
!
!////////////////////////////////////////////////////////////////////////
!
!>The FTKeyObjectPairClass is for use by the FTDictionary Class and will
!>generally not be interacted with by the user.
!>
      Module FTKeyObjectPairClass
         USE FTObjectClass
         IMPLICIT NONE  
!
!    -----------------
!    Module constants:
!    -----------------
!
         INTEGER, PARAMETER, PUBLIC  :: FTDICT_KWD_STRING_LENGTH    = 64
!
!        ----------
!        Class type
!        ----------
!
         TYPE, EXTENDS(FTObject) :: FTKeyObjectPair
            CHARACTER(LEN=FTDICT_KWD_STRING_LENGTH) :: keyString
            CLASS(FTObject) , POINTER               :: valueObject => NULL()
!
!           --------            
            CONTAINS
!           --------
!          
            PROCEDURE :: initWithObjectAndKey
            PROCEDURE :: destruct         => destructFTKeyObjectPair
            PROCEDURE :: description      => FTKeyObjectPairDescription
            PROCEDURE :: printDescription => printFTKeyObjectPairDescription
            
            PROCEDURE :: key
            PROCEDURE :: object

         END TYPE FTKeyObjectPair
      
      INTERFACE release
         MODULE PROCEDURE releaseFTKeyObjectPair
      END INTERFACE  
!
!        ========       
         CONTAINS  
!        ========       
!
!@mark -
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE initWithObjectAndKey(self,v,key)
            IMPLICIT NONE
            CLASS(FTKeyObjectPair)      :: self
            CHARACTER(LEN=*)            :: key
            CLASS(FTObject) , POINTER   :: v
            
            CALL self % FTObject % init
            
            self % keyString   = key
            self % valueObject => v
            
            CALL self % valueObject % retain
            
         END SUBROUTINE initWithObjectAndKey
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE destructFTKeyObjectPair(self)
            IMPLICIT NONE
            CLASS(FTKeyObjectPair) :: self
            
            self % keyString = ""
            CALL releaseFTObject(self % valueObject)
!
!           ----------------------------------------------------
!           Call superclass destructor after processing subclass
!           specific items
!           ----------------------------------------------------
!
            CALL self % FTObject % destruct()
             
         END SUBROUTINE destructFTKeyObjectPair
!
!------------------------------------------------
!> Public, generic name: release(self)
!>
!> Call release(self) on an object to release control
!> of an object. If its reference count is zero, then 
!> it is deallocated.
!------------------------------------------------
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE releaseFTKeyObjectPair(self)  
         IMPLICIT NONE
         CLASS(FTKeyObjectPair) , POINTER :: self
         CLASS(FTObject)        , POINTER :: obj
         obj => self
         CALL releaseFTObject(self = obj)
         IF ( .NOT. ASSOCIATED(obj) )     THEN
            self => NULL() 
         END IF      
      END SUBROUTINE releaseFTKeyObjectPair
!
!//////////////////////////////////////////////////////////////////////// 
! 
         CHARACTER(LEN=DESCRIPTION_CHARACTER_LENGTH) FUNCTION FTKeyObjectPairDescription(self)  
            IMPLICIT NONE  
            CLASS(FTKeyObjectPair) :: self
            
            WRITE(FTKeyObjectPairDescription,*)"(" , TRIM(self % keyString)  , "," &
                                           , TRIM(self % valueObject % description()) , ")"
             
         END FUNCTION FTKeyObjectPairDescription
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE printFTKeyObjectPairDescription(self,iUnit)  
            IMPLICIT NONE  
            CLASS(FTKeyObjectPair) :: self
            INTEGER                :: iUnit
            
            WRITE(iUnit,*) "{"

            WRITE(iUnit,'(6x,A,A3)',ADVANCE = "NO") TRIM(self % keyString)  , " = " 
            CALL self % valueObject % printDescription(iUnit)
            
            WRITE(iUnit,*) "}"
             
         END SUBROUTINE printFTKeyObjectPairDescription    
!
!//////////////////////////////////////////////////////////////////////// 
! 
         CHARACTER(LEN=FTDICT_KWD_STRING_LENGTH) FUNCTION key(self)  
            IMPLICIT NONE  
            CLASS(FTKeyObjectPair) :: self
            key = self % keyString
         END FUNCTION key         
!
!//////////////////////////////////////////////////////////////////////// 
! 
         FUNCTION object(self)  
            IMPLICIT NONE  
            CLASS(FTKeyObjectPair)   :: self
            CLASS(FTObject), POINTER :: object
            
            object => self % valueObject
          
         END FUNCTION object    
          
      END Module FTKeyObjectPairClass  
!
!////////////////////////////////////////////////////////////////////////
!
!@mark -
!>
!>A dictionary is a special case of a hash table that stores key-value pairs. 
!>
!>It is an
!>example of what is called an ``associative container''. In the implementation of FTObjectLibrary,
!>the value can be any subclass of FTObject and the key is a character variable. The library
!>includes the base dictionary that can store and retrieve any subclass of FTObject. It also includes a
!>subclass that is designed to store and retrieve FTValue objects.
!> 
!>FTDictionary (Inherits from FTObject)
!>
!>###Initialization
!>
!>         CLASS(FTDictionary), POINTER :: dict
!>         ALLOCATE(dict)
!>         CALL dict % initWithSize(N) ! N = size of dictionary. Should be power of two
!>
!>###Adding entries
!>
!>         CLASS(FTDictionary), POINTER :: dict
!>         CLASS(FTObject)    , POINTER :: obj
!>         CHARACTER(LEN=N)             :: key
!>         obj => r                            ! r is subclass of FTObject
!>         CALL dict % addObjectForKey(obj,key)
!>
!>###Accessing entries
!>
!>         obj => dict % objectForKey(key)
!>         CALL cast(obj,v) ! v is the type of object to be extracted
!>
!>###Destruction
!>   
!>         CALL release(dict) ! Pointer
!>         call dict % destruct() ! Stack variable
!>###Accessing an object
!>
!>           TYPE(FTObject) :: obj
!>           obj => dict % objectForKey(key)
!>
!>###Converting a base class pointer to a dictionary
!>           dict =>  dictionaryFromObject(obj)
!>
!>###Getting all of the keys
!>           CHARACTER(LEN=FTDICT_KWD_STRING_LENGTH), POINTER :: keys(:)
!>           keys =>  dict % allKeys()
!>(The target of the pointer must be deallocated by the caller)
!>###Getting all of the objects
!>           CLASS(FTMutableObjectArray), POINTER :: objectArray
!>           objectArray =>  dict % allObjects() ! The array is owned by the caller.
!>(The target of the pointer must be released by the caller)
!>
      Module FTDictionaryClass
         USE FTKeyObjectPairClass
         USE FTLinkedListClass
         USE FTLinkedListIteratorClass
         USE FTMutableObjectArrayClass
         IMPLICIT NONE  
         
         TYPE, EXTENDS(FTObject) :: FTDictionary
            INTEGER                                   :: numberOfEntries
            LOGICAL                                   :: isCaseSensitive
            TYPE(FTLinkedList), DIMENSION(:), POINTER :: entries => NULL()
!
!           --------
            CONTAINS
!           --------
!
            PROCEDURE :: initWithSize
            PROCEDURE :: init
            PROCEDURE :: setCaseSensitive
            PROCEDURE :: caseSensitive
            PROCEDURE :: allKeys
            PROCEDURE :: allObjects
            PROCEDURE :: destruct => destructFTDictionary
            PROCEDURE :: addObjectForKey
            procedure :: removeObjectForKey
            PROCEDURE :: description => FTDictionaryDescription
            PROCEDURE :: printDescription => printFTDictionaryDescription
            PROCEDURE :: objectForKey
            PROCEDURE :: containsKey
            PROCEDURE :: className => dictionaryClassName
            PROCEDURE :: COUNT
         END TYPE FTDictionary
         
         INTERFACE cast
            MODULE PROCEDURE castToDictionary
         END INTERFACE cast
         
         INTERFACE release
            MODULE PROCEDURE releaseFTDictionary 
         END INTERFACE  
!
!        ========         
         CONTAINS  
!        ========  
!       
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE init(self)  
            IMPLICIT NONE  
            CLASS(FTDictionary) :: self
            CALL initWithSize(self,16)
         END SUBROUTINE init
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE initWithSize(self,sze)  
!
!           ----------------------
!           Designated initializer
!           ----------------------
!
            IMPLICIT NONE
            CLASS(FTDictionary) :: self
            INTEGER, INTENT(in) :: sze
            INTEGER             :: i
            
            CALL self % FTObject % init()
            
            self % isCaseSensitive = .true.
            self % numberOfEntries = 0
!
!           --------------------------------
!           Create the array of linked lists
!           --------------------------------
!
            ALLOCATE(self % entries(sze))
            DO i = 1, sze
               CALL self % entries(i) % init()
            END DO    
            
         END SUBROUTINE initWithSize
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE destructFTDictionary(self)  
            IMPLICIT NONE  
            CLASS(FTDictionary) :: self
           
            INTEGER :: i
            
            DO i = 1, SIZE(self % entries)
               CALL self % entries(i) % destruct()
            END DO

            DEALLOCATE(self % entries)
            self % entries => NULL()
!
!           ----------------------------------------------------
!           Call superclass destructor after processing subclass
!           specific items
!           ----------------------------------------------------
!
           CALL self % FTObject % destruct()
            
         END SUBROUTINE destructFTDictionary    
!
!------------------------------------------------
!> Public, generic name: release(self)
!>
!> Call release(self) on an object to release control
!> of an object. If its reference count is zero, then 
!> it is deallocated.
!------------------------------------------------
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE releaseFTDictionary(self)  
         IMPLICIT NONE
         TYPE(FTDictionary) , POINTER :: self
         CLASS(FTObject)    , POINTER :: obj
         
         IF(.NOT. ASSOCIATED(self)) RETURN
         
         obj => self
         CALL releaseFTObject(self = obj)
         IF ( .NOT. ASSOCIATED(obj) )     THEN
            self => NULL() 
         END IF      
      END SUBROUTINE releaseFTDictionary
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE setCaseSensitive(self,bool)  
            IMPLICIT NONE  
            CLASS(FTDictionary) :: self
            LOGICAL             :: bool
            self % isCaseSensitive = bool
         END SUBROUTINE setCaseSensitive  
!
!//////////////////////////////////////////////////////////////////////// 
! 
         LOGICAL FUNCTION caseSensitive(self)  
            IMPLICIT NONE  
            CLASS(FTDictionary) :: self
            caseSensitive = self % isCaseSensitive
         END FUNCTION caseSensitive      
!
!//////////////////////////////////////////////////////////////////////// 
! 
         INTEGER FUNCTION COUNT(self)  
            IMPLICIT NONE  
            CLASS(FTDictionary) :: self
            COUNT = self % numberOfEntries
         END FUNCTION COUNT    
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE addObjectForKey(self,object,key)
            IMPLICIT NONE  
            CLASS(FTDictionary)                 :: self
            CLASS(FTObject)       , POINTER     :: object
            CHARACTER(LEN=*)                    :: key
            CLASS(FTKeyObjectPair), POINTER     :: pair => NULL()
            CLASS(FTObject)       , POINTER     :: ptr => NULL()
            INTEGER                             :: h
            
            INTEGER, EXTERNAL                   :: b3hs_hash_key_jenkins
            
            h = b3hs_hash_key_jenkins(key,SIZE(self % entries))
            
            ALLOCATE(pair)
            CALL pair % initWithObjectAndKey(object,key)
            ptr => pair
            CALL self % entries(h) % add(ptr)
            self % numberOfEntries = self % numberOfEntries + 1
            
         END SUBROUTINE addObjectForKey

         SUBROUTINE removeObjectForKey(self,key)
            IMPLICIT NONE  
            CLASS(FTDictionary)                 :: self
            CHARACTER(LEN=*)                    :: key
!
!           ---------------
!           Local variables
!           ---------------
!
            CLASS(FTObject)       , POINTER     :: obj => NULL()
            INTEGER                             :: h
            CLASS(FTLinkedListRecord)     , POINTER :: listRecordPtr => NULL()
            CHARACTER(LEN=FTDICT_KWD_STRING_LENGTH) :: keyString
            
            INTEGER, EXTERNAL                   :: b3hs_hash_key_jenkins
            
            h = b3hs_hash_key_jenkins(key,SIZE(self % entries))

            obj => NULL()
            
            listRecordPtr => self % entries(h) % head
            DO WHILE(ASSOCIATED(listRecordPtr))
!
!              --------------------------------------------
!              The self % entries(h)'s recordObject is a FTKeyObjectPair
!              --------------------------------------------
!
               SELECT TYPE (pair => listRecordPtr % recordObject)
                  TYPE is (FTKeyObjectPair)
                     keyString = pair % key()
                     IF ( TRIM(keyString) == TRIM(key) )     THEN
                        obj => pair
                        EXIT 
                     END IF 
                  CLASS DEFAULT
               END SELECT
               listRecordPtr  => listRecordPtr % next
            END DO    

            CALL self % entries(h) % remove(obj)
            self % numberOfEntries = self % numberOfEntries - 1
            
         END SUBROUTINE removeObjectForKey
!
!//////////////////////////////////////////////////////////////////////// 
! 
         FUNCTION objectForKey(self,key)  
            IMPLICIT NONE  
            CLASS(FTDictionary)                  :: self
            CHARACTER(LEN=*)                     :: key
            CLASS(FTObject)       , POINTER      :: objectForKey
            INTEGER                              :: h
            
            INTEGER, EXTERNAL                    :: b3hs_hash_key_jenkins
           
            objectForKey => NULL()
            IF(self % COUNT() == 0)     RETURN 
            
            h = b3hs_hash_key_jenkins(key,SIZE(self % entries))
            IF ( self % entries(h) % COUNT() > 0 )     THEN
               objectForKey => objectForKeyInList(key,self % entries(h))
            END IF 
 
         END FUNCTION ObjectForKey
!
!//////////////////////////////////////////////////////////////////////// 
! 
         FUNCTION containsKey(self,key)  RESULT(r)
            IMPLICIT NONE  
            CLASS(FTDictionary)      :: self
            CHARACTER(LEN=*)         :: key
            CLASS(FTObject), POINTER :: obj
            LOGICAL                  :: r
           
            IF ( ASSOCIATED( self % objectForKey(key)) )     THEN
               r = .TRUE. 
            ELSE 
               r = .FALSE. 
            END IF 
 
         END FUNCTION containsKey
!
!//////////////////////////////////////////////////////////////////////// 
! 
         FUNCTION objectForKeyInList(key,list)  
            IMPLICIT NONE  
            CHARACTER(LEN=*)            :: key
            CLASS(FTLinkedList)         :: list
            CLASS(FTObject), POINTER    :: objectForKeyInList
            
            CLASS(FTLinkedListRecord)     , POINTER :: listRecordPtr => NULL()
            CHARACTER(LEN=FTDICT_KWD_STRING_LENGTH) :: keyString

            objectForKeyInList => NULL()
            
            listRecordPtr => list % head
            DO WHILE(ASSOCIATED(listRecordPtr))
!
!              --------------------------------------------
!              The list's recordObject is a FTKeyObjectPair
!              --------------------------------------------
!
               SELECT TYPE (pair => listRecordPtr % recordObject)
                  TYPE is (FTKeyObjectPair)
                     keyString = pair % key()
                     IF ( TRIM(keyString) == TRIM(key) )     THEN
                        objectForKeyInList => pair % object()
                        EXIT 
                     END IF 
                  CLASS DEFAULT
               END SELECT
               listRecordPtr  => listRecordPtr % next
            END DO    

         END FUNCTION objectForKeyInList
!
!//////////////////////////////////////////////////////////////////////// 
! 
         CHARACTER(LEN=DESCRIPTION_CHARACTER_LENGTH) FUNCTION FTDictionaryDescription(self)  
            IMPLICIT NONE  
            CLASS(FTDictionary) :: self
            CHARACTER(LEN=DESCRIPTION_CHARACTER_LENGTH) :: s
            
            INTEGER :: i
            FTDictionaryDescription = ""

            IF(SELF % COUNT() == 0) RETURN
            
            DO i = 1, SIZE(self % entries)
               s = self % entries(i) % description()
               IF ( LEN_TRIM(s) > 0 )     THEN
                  FTDictionaryDescription =  TRIM(FTDictionaryDescription) // &
                                             TRIM(self % entries(i) % description()) // &
                                             CHAR(13)
               END IF 
            END DO
            
         END FUNCTION FTDictionaryDescription    
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE printFTDictionaryDescription(self,iUnit)  
            IMPLICIT NONE  
            CLASS(FTDictionary) :: self
            INTEGER             :: iUnit
            
            INTEGER :: i

            IF(SELF % COUNT() == 0) THEN
               WRITE(iUnit,*) "Empty Dictionary"
            END IF  
            
            DO i = 1, SIZE(self % entries)
               CALL self % entries(i) % printDescription(iUnit)
            END DO
            
         END SUBROUTINE printFTDictionaryDescription
!
!//////////////////////////////////////////////////////////////////////// 
! 
         FUNCTION AllObjects(self) RESULT(objectArray)
            IMPLICIT NONE  
!
!           ---------
!           Arguments
!           ---------
!
            CLASS(FTDictionary)                  :: self
            CLASS(FTMutableObjectArray), POINTER :: objectArray
!
!           ---------------
!           Local Variables
!           ---------------
!
            INTEGER                                 :: i
            CLASS(FTLinkedListRecord)     , POINTER :: listRecordPtr => NULL()
            CLASS(FTObject)               , POINTER :: obj           => NULL()
            CHARACTER(LEN=FTDICT_KWD_STRING_LENGTH) :: keyString
!
!           --------------------------------------------
!           Allocate a pointer to the object array to be
!           returned with refCount = 1
!           --------------------------------------------
!
            ALLOCATE(objectArray)
            CALL objectArray % initWithSize(arraySize = self % COUNT())
            
            DO i = 1, SIZE(self % entries)
               listRecordPtr => self % entries(i) % head
               DO WHILE(ASSOCIATED(listRecordPtr))
!
!                 --------------------------------------------
!                 The list's recordObject is a FTKeyObjectPair
!                 --------------------------------------------
!
                  SELECT TYPE (pair => listRecordPtr % recordObject)
                     TYPE is (FTKeyObjectPair)
                        keyString = pair % key()
                        obj  => pair % object()
                        CALL objectArray % addObject(obj)
                     CLASS DEFAULT
                  END SELECT
                  listRecordPtr  => listRecordPtr % next
               END DO    
            END DO  
            
         END FUNCTION AllObjects
!
!//////////////////////////////////////////////////////////////////////// 
! 
         FUNCTION AllKeys(self) RESULT(keys)
            IMPLICIT NONE  
!
!           ---------
!           Arguments
!           ---------
!
            CLASS(FTDictionary)                              :: self
            CHARACTER(LEN=FTDICT_KWD_STRING_LENGTH), POINTER :: keys(:)
!
!           ---------------
!           Local Variables
!           ---------------
!
            INTEGER                                 :: i, c
            CLASS(FTLinkedListRecord)     , POINTER :: listRecordPtr => NULL()
            CHARACTER(LEN=FTDICT_KWD_STRING_LENGTH) :: keyString
!
!           --------------------------------------------
!           Allocate a pointer to the object array to be
!           returned with refCount = 1
!           --------------------------------------------
!
            ALLOCATE(keys(self % COUNT()))
            
            c = 1
            DO i = 1, SIZE(self % entries)
               listRecordPtr => self % entries(i) % head
               DO WHILE(ASSOCIATED(listRecordPtr))
!
!                 --------------------------------------------
!                 The list's recordObject is a FTKeyObjectPair
!                 --------------------------------------------
!
                  SELECT TYPE (pair => listRecordPtr % recordObject)
                     TYPE is (FTKeyObjectPair)
                        keyString = pair % key()
                        keys(c)   = keyString
                     CLASS DEFAULT
                  END SELECT
                  c = c + 1
                  listRecordPtr  => listRecordPtr % next
               END DO    
            END DO  
            
         END FUNCTION AllKeys
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE castToDictionary(obj,cast) 
!
!     -----------------------------------------------------
!     Cast the base class FTObject to the FTException class
!     -----------------------------------------------------
!
         IMPLICIT NONE  
         CLASS(FTObject)    , POINTER :: obj
         CLASS(FTDictionary), POINTER :: cast
         
         cast => NULL()
         SELECT TYPE (e => obj)
            TYPE is (FTDictionary)
               cast => e
            CLASS DEFAULT
               
         END SELECT
         
      END SUBROUTINE castToDictionary
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION dictionaryFromObject(obj) RESULT(cast)
!
!     -----------------------------------------------------
!     Cast the base class FTObject to the FTException class
!     -----------------------------------------------------
!
         IMPLICIT NONE  
         CLASS(FTObject)    , POINTER :: obj
         CLASS(FTDictionary), POINTER :: cast
         
         cast => NULL()
         SELECT TYPE (e => obj)
            TYPE is (FTDictionary)
               cast => e
            CLASS DEFAULT
               
         END SELECT
         
      END FUNCTION dictionaryFromObject
!
!//////////////////////////////////////////////////////////////////////// 
! 
!      -----------------------------------------------------------------
!> Class name returns a string with the name of the type of the object
!>
!>  ### Usage:
!>
!>        PRINT *,  obj % className()
!>        if( obj % className = "FTDictionary")
!>
      FUNCTION dictionaryClassName(self)  RESULT(s)
         IMPLICIT NONE  
         CLASS(FTDictionary)                        :: self
         CHARACTER(LEN=CLASS_NAME_CHARACTER_LENGTH) :: s
         
         s = "FTDictionary"
 
      END FUNCTION dictionaryClassName

      END Module FTDictionaryClass    
      