!
!////////////////////////////////////////////////////////////////////////
!
!      FTLinkedListClass.f90
!      Created: January 7, 2013 2:56 PM 
!      By: David Kopriva  
!
!
!
!////////////////////////////////////////////////////////////////////////
!
!@mark -
!
!>FTLinkedListRecord is the record type (object and next) for the
!>LinkedList class.
!>
!>One will generally not instantiate a record oneself. They are 
!>created automatically when one adds an object to a linked list.
!>
      Module FTLinkedListRecordClass 
      USE FTObjectClass
      IMPLICIT NONE 
!
!     -----------------------------
!     Record class for linked lists
!     -----------------------------
!
      TYPE, EXTENDS(FTObject) :: FTLinkedListRecord
      
         CLASS(FTObject)          , POINTER :: recordObject => NULL()
         CLASS(FTLinkedListRecord), POINTER :: next => NULL(), previous => NULL()
!
!        ========         
         CONTAINS
!        ========
!
         PROCEDURE :: initWithObject
         PROCEDURE :: destruct         => destructFTLinkedListRecord
         PROCEDURE :: printDescription => printFTLinkedRecordDescription
         PROCEDURE :: className        => llRecordClassName
         
      END TYPE FTLinkedListRecord
      
      INTERFACE release
         MODULE PROCEDURE releaseFTLinkedListRecord
      END INTERFACE  
!
!     ----------
!     Procedures
!     ----------
!
      CONTAINS 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE initWithObject(self,obj) 
         IMPLICIT NONE 
         CLASS(FTLinkedListRecord) :: self
         CLASS(FTObject), POINTER  :: obj
!
!        -------------------------------
!        Always call the superclass init
!        -------------------------------
!
         CALL self % FTObject % init()
!
!        ------------------------
!        Subclass initializations
!        ------------------------
!
         CALL obj % retain()
         
         self % recordObject => obj
         self % next         => NULL()
         self % previous     => NULL()
         
      END SUBROUTINE initWithObject
!
!////////////////////////////////////////////////////////////////////////
!
!< The destructor must only be called from within subclass destructors
!
      SUBROUTINE destructFTLinkedListRecord(self) 
         IMPLICIT NONE
         CLASS(FTLinkedListRecord) :: self
         
         IF ( ASSOCIATED(self % recordObject) ) CALL releaseFTObject(self % recordObject)
         self % next     => NULL()
         self % previous => NULL()
!
!        ------------------------------------------
!        Always call the superclass destructor here
!        at the end of the subclass destructor.
!        ------------------------------------------
!
         CALL self % FTObject % destruct()
        
      END SUBROUTINE destructFTLinkedListRecord
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE releaseFTLinkedListRecord(self)  
         IMPLICIT NONE
         CLASS(FTLinkedListRecord) , POINTER :: self
         CLASS(FTObject)           , POINTER :: obj
         
         IF(.NOT. ASSOCIATED(self)) RETURN
         
         obj => self
         CALL releaseFTObject(self = obj)
         IF ( .NOT. ASSOCIATED(obj) )     THEN
            self => NULL() 
         END IF      
      END SUBROUTINE releaseFTLinkedListRecord
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE printFTLinkedRecordDescription(self,iUnit)  
         IMPLICIT NONE  
         CLASS(FTLinkedListRecord) :: self
         INTEGER                   :: iUnit

         IF ( ASSOCIATED(self % recordObject) )     THEN
            CALL self % recordObject % printDescription(iUnit)
         END IF 
         
      END SUBROUTINE printFTLinkedRecordDescription
!
!//////////////////////////////////////////////////////////////////////// 
! 
!      -----------------------------------------------------------------
!> Class name returns a string with the name of the type of the object
!>
!>  ### Usage:
!>
!>        PRINT *,  obj % className()
!>        if( obj % className = "FTLinkedListRecord")
!>
      FUNCTION llRecordClassName(self)  RESULT(s)
         IMPLICIT NONE  
         CLASS(FTLinkedListRecord)                  :: self
         CHARACTER(LEN=CLASS_NAME_CHARACTER_LENGTH) :: s
         
         s = "FTLinkedListRecord"
 
      END FUNCTION llRecordClassName

      
      END MODULE FTLinkedListRecordClass  
!@mark -
!
!
!     --------------------------------------------------
!     Implements the basics of a linked list of objects
!     --------------------------------------------------
!
!>
!>FTLinkedList is a container class that stores objects in a linked list.
!>
!>Inherits from FTObjectClass
!>
!>##Definition (Subclass of FTObject):
!>
!>         TYPE(FTLinkedList) :: list
!>
!>#Usage:
!>
!>##Initialization
!>
!>         CLASS(FTLinkedList), POINTER :: list
!>         ALLOCATE(list)
!>         CALL list % init
!>
!>##Adding objects
!>
!>         CLASS(FTLinkedList), POINTER :: list, listToAdd
!>         CLASS(FTObject)    , POINTER :: objectPtr
!>
!>         objectPtr => r                ! r is subclass of FTObject
!>         CALL list % Add(objectPtr)    ! Pointer is retained by list
!>         CALL release(r)               ! If caller relinquishes ownership
!>
!>         CALL list % addObjectsFromList(listToAdd)
!>
!>##Inserting objects
!>
!>         CLASS(FTLinkedList)      , POINTER :: list
!>         CLASS(FTObject)          , POINTER :: objectPtr, obj
!>         CLASS(FTLinkedListRecord), POINTER :: record
!>
!>         objectPtr => r                                        ! r is subclass of FTObject
!>         CALL list % insertObjectAfterRecord(objectPtr,record) ! Pointer is retained by list
!>         CALL release(r)                                       ! If caller reliquishes ownership
!>
!>         objectPtr => r                                     ! r is subclass of FTObject
!>         CALL list % insertObjectAfterObject(objectPtr,obj) ! Pointer is retained by list
!>         CALL release(r)                                    ! If caller reliquishes ownership
!>
!>##Removing objects
!>
!>         CLASS(FTLinkedList), POINTER :: list
!>         CLASS(FTObject)    , POINTER :: objectPtr
!>         objectPtr => r                 ! r is subclass of FTObject
!>         CALL list % remove(objectPtr)
!>
!>##Getting all objects as an object array
!>
!>         CLASS(FTLinkedList)        , POINTER :: list
!>         CLASS(FTMutableObjectArray), POINTER :: array
!>         array => list % allObjects() ! Array has refCount = 1
!>
!>##Counting the number of objects in the list
!>
!>         n = list % count()
!>
!>##Destruction
!>   
!>         CALL release(list) [Pointers]
!>         CALL list % destruct() [Non Pointers]
!>!
      Module FTLinkedListClass
!      
      USE FTLinkedListRecordClass
      USE FTMutableObjectArrayClass
      IMPLICIT NONE 
!
!     -----------------
!     Class object type
!     -----------------
!
      TYPE, EXTENDS(FTObject) :: FTLinkedList
      
         CLASS(FTLinkedListRecord), POINTER :: head => NULL(), tail => NULL()
         INTEGER                            :: nRecords
         LOGICAL                            :: isCircular_
!
!        ========         
         CONTAINS
!        ========
!
         PROCEDURE :: init             => initFTLinkedList
         PROCEDURE :: add              
         PROCEDURE :: remove           => removeObject
         PROCEDURE :: reverse          => reverseLinkedList
         PROCEDURE :: removeRecord     => removeLinkedListRecord
         PROCEDURE :: destruct         => destructFTLinkedList
         PROCEDURE :: count            => numberOfRecords
         PROCEDURE :: description      => FTLinkedListDescription
         PROCEDURE :: printDescription => printFTLinkedListDescription
         PROCEDURE :: className        => linkedListClassName
         PROCEDURE :: allObjects       => allLinkedListObjects
         PROCEDURE :: removeAllObjects => removeAllLinkedListObjects
         PROCEDURE :: addObjectsFromList
         PROCEDURE :: makeCircular
         PROCEDURE :: isCircular
         PROCEDURE :: insertObjectAfterRecord
         PROCEDURE :: insertObjectAfterObject
         
      END TYPE FTLinkedList
      
      INTERFACE cast
         MODULE PROCEDURE castObjectToLinkedList
      END INTERFACE cast
      
      INTERFACE release
         MODULE PROCEDURE releaseFTLinkedList 
      END INTERFACE  
      
!
!     ----------
!     Procedures
!     ----------
!
      CONTAINS 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE initFTLinkedList(self) 
         IMPLICIT NONE 
         CLASS(FTLinkedList) :: self
!
!        -------------------------------
!        Always call the superclass init
!        -------------------------------
!
         CALL self % FTObject % init()
!
!        --------------------------------------
!        Then call the subclass initializations
!        --------------------------------------
!
         self % nRecords    = 0
         self % isCircular_ = .FALSE.
         
         self % head => NULL(); self % tail => NULL()
         
      END SUBROUTINE initFTLinkedList
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE add(self,obj)
         IMPLICIT NONE 
         CLASS(FTLinkedList)                :: self
         CLASS(FTObject)          , POINTER :: obj
         CLASS(FTLinkedListRecord), POINTER :: newRecord => NULL()
         
         ALLOCATE(newRecord)
         CALL newRecord % initWithObject(obj)
         
         IF ( .NOT.ASSOCIATED(self % head) )     THEN
            self % head => newRecord
         ELSE
            self % tail % next   => newRecord
            newRecord % previous => self % tail
         END IF
         
         self % tail => newRecord
         self % nRecords = self % nRecords + 1
         
      END SUBROUTINE add 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE addObjectsFromList(self,list)
         IMPLICIT NONE 
         CLASS(FTLinkedList)                :: self
         CLASS(FTLinkedList)      , POINTER :: list
         CLASS(FTLinkedListRecord), POINTER :: recordPtr => NULL()
         CLASS(FtObject)          , POINTER :: obj       => NULL()
         LOGICAL                            :: circular
         
         IF(.NOT.ASSOCIATED(list % head)) RETURN
         
         circular = list % isCircular()
         CALL list % makeCircular(.FALSE.)

         recordPtr => list % head
         DO WHILE(ASSOCIATED( recordPtr ))
            obj => recordPtr % recordObject
            CALL self % add(obj)
            
            recordPtr => recordPtr % next
         END DO
         
         CALL list % makeCircular(circular) 
         
      END SUBROUTINE addObjectsFromList 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE insertObjectAfterRecord(self,obj,after)
         IMPLICIT NONE 
         CLASS(FTLinkedList)                :: self
         CLASS(FTObject)          , POINTER :: obj
         CLASS(FTLinkedListRecord), POINTER :: newRecord   => NULL()
         CLASS(FTLinkedListRecord), POINTER :: after, next => NULL()
         
         ALLOCATE(newRecord)
         CALL newRecord % initWithObject(obj)
         
         next                 => after % next
         newRecord % next     => next
         newRecord % previous => after
         after % next         => newRecord
         next % previous      => newRecord
         
         IF ( .NOT.ASSOCIATED( newRecord % next ) )     THEN
            self % tail => newRecord 
         END IF 
         
         self % nRecords = self % nRecords + 1
         
      END SUBROUTINE insertObjectAfterRecord 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE insertObjectAfterObject(self,obj,after)
         IMPLICIT NONE 
         CLASS(FTLinkedList)                :: self
         CLASS(FTObject)          , POINTER :: obj, after
         
         CLASS(FTLinkedListRecord), POINTER :: current => NULL(), previous => NULL()
                  
         IF ( .NOT.ASSOCIATED(self % head) )     THEN
            CALL self % add(obj)
            RETURN 
         END IF 
         
         current  => self % head
         previous => NULL()
!
!        -------------------------------------------------------------
!        Find the object in the list by a linear search and 
!        add the new object after it.
!        It will be deallocated if necessary.
!        -------------------------------------------------------------
!
         DO WHILE (ASSOCIATED(current))
         
            IF ( ASSOCIATED(current % recordObject,after) )     THEN
               CALL self % insertObjectAfterRecord(obj = obj,after = current)
               RETURN 
            END IF 
            
            previous => current
            current  => current % next
         END DO
         
      END SUBROUTINE insertObjectAfterObject 
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE makeCircular(self,circular)  
         IMPLICIT NONE  
         CLASS(FTLinkedList) :: self
         LOGICAL             :: circular
         
         IF ( circular )     THEN
            self % head % previous => self % tail
            self % tail % next     => self % head
            self % isCircular_ = .TRUE.
         ELSE
            self % head % previous => NULL()
            self % tail % next     => NULL()
            self % isCircular_ = .FALSE.
         END IF 
      END SUBROUTINE makeCircular
!
!//////////////////////////////////////////////////////////////////////// 
! 
      LOGICAL FUNCTION isCircular(self)  
         IMPLICIT NONE  
         CLASS(FTLinkedList) :: self
         isCircular = self % isCircular_
      END FUNCTION isCircular
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE removeObject1(self,obj)
         IMPLICIT NONE 
         CLASS(FTLinkedList)                :: self
         CLASS(FTObject)          , POINTER :: obj
         
         CLASS(FTLinkedListRecord), POINTER :: current => NULL(), previous => NULL()
                  
         IF ( .NOT.ASSOCIATED(self % head) )     RETURN
         
         current  => self % head
         previous => NULL()
!
!        -------------------------------------------------------------
!        Find the object in the list by a linear search and remove it.
!        It will be deallocated if necessary.
!        -------------------------------------------------------------
!
         DO WHILE (ASSOCIATED(current))
            IF ( ASSOCIATED(current % recordObject,obj) )     THEN
            
               IF ( ASSOCIATED(previous) )     THEN
                  previous % next => current % next
               ELSE
                  self % head => current % next
               END IF 
               
               IF ( ASSOCIATED(current,self % tail) )     THEN
                  self % tail => previous 
               END IF 
               
               CALL releaseFTLinkedListRecord(current)
               
               self % nRecords = self % nRecords - 1
               EXIT
            END IF 
            
            previous => current
            current  => current % next
         END DO
         
      END SUBROUTINE removeObject1 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE removeObject(self,obj)
         IMPLICIT NONE 
         CLASS(FTLinkedList)                :: self
         CLASS(FTObject)          , POINTER :: obj
         
         CLASS(FTLinkedListRecord), POINTER :: current => NULL(), previous => NULL()
                  
         IF ( .NOT.ASSOCIATED(self % head) )     RETURN
         
         current  => self % head
         previous => NULL()
!
!        -------------------------------------------------------------
!        Find the object in the list by a linear search and remove it.
!        It will be deallocated if necessary.
!        -------------------------------------------------------------
!
         DO WHILE (ASSOCIATED(current))
         
            IF ( ASSOCIATED(current % recordObject,obj) )     THEN
               CALL self % removeRecord(current)
               EXIT
            END IF 
            
            previous => current
            current  => current % next
         END DO
         
      END SUBROUTINE removeObject 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE removeLinkedListRecord(self,listRecord)
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(FTLinkedList)                :: self
         CLASS(FTLinkedListRecord), POINTER :: listRecord
!
!        ---------------
!        Local variables
!        ---------------
!
         CLASS(FTLinkedListRecord), POINTER :: previous => NULL(), next => NULL()
!
!        ---------------------------------------------------
!        Turn cirularity off and then back on
!        to work around an what appears to be an
!        ifort bug testing the association of two pointers. 
!        ---------------------------------------------------
!
         LOGICAL :: circ
         circ = self % isCircular()
         IF(circ) CALL self % makeCircular(.FALSE.)
         
         previous => listRecord % previous
         next     => listRecord % next
         
         IF ( .NOT.ASSOCIATED(listRecord % previous) )     THEN
            self % head => next
            IF ( ASSOCIATED(next) )     THEN
               self % head % previous => NULL() 
            END IF  
         END IF 
         
         IF ( .NOT.ASSOCIATED(listRecord % next) )     THEN
            self % tail => previous
            IF ( ASSOCIATED(previous) )     THEN
               self % tail % next => NULL() 
            END IF
         END IF 
         
         IF ( ASSOCIATED(previous) .AND. ASSOCIATED(next) )     THEN
            previous % next => next
            next % previous => previous 
         END IF 
         
         CALL releaseFTLinkedListRecord(listRecord)
         
         self % nRecords = self % nRecords - 1
         IF(circ) CALL self % makeCircular(.TRUE.)
         
      END SUBROUTINE removeLinkedListRecord
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE removeAllLinkedListObjects(self)  
         IMPLICIT NONE
         CLASS(FTLinkedList)                :: self
         CLASS(FTLinkedListRecord), POINTER :: listRecord => NULL(), tmp => NULL()
         LOGICAL                            :: circular

         IF(.NOT.ASSOCIATED(self % head)) RETURN 
         
         circular = self % isCircular()
         CALL self % makeCircular(.FALSE.)
         
         listRecord => self % head
         DO WHILE (ASSOCIATED(listRecord))

            tmp => listRecord % next

            CALL releaseFTLinkedListRecord(listRecord)

            IF(.NOT. ASSOCIATED(listRecord)) THEN
               self % nRecords = self % nRecords - 1
            END IF
            listRecord => tmp
         END DO

         self % head => NULL(); self % tail => NULL()
         
      END SUBROUTINE removeAllLinkedListObjects
!
!//////////////////////////////////////////////////////////////////////// 
! 
      INTEGER FUNCTION numberOfRecords(self)  
         IMPLICIT NONE  
         CLASS(FTLinkedList) :: self
          numberOfRecords = self % nRecords
      END FUNCTION numberOfRecords     
!
!////////////////////////////////////////////////////////////////////////
!
!< The destructor must only be called from within the destructors of subclasses
!> It is automatically called by release().
!
      SUBROUTINE destructFTLinkedList(self) 
         IMPLICIT NONE
         CLASS(FTLinkedList)                :: self

         CALL self % removeAllObjects()
!
!        ------------------------------------------
!        Always call the superclass destructor here
!        at the end of the subclass destructor.
!        ------------------------------------------
!
         CALL self % FTObject % destruct()

      END SUBROUTINE destructFTLinkedList
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
      SUBROUTINE releaseFTLinkedList(self)  
         IMPLICIT NONE
         TYPE (FTLinkedList), POINTER :: self
         CLASS(FTObject)    , POINTER :: obj
         
         IF(.NOT. ASSOCIATED(self)) RETURN

         obj => self
         CALL releaseFTObject(self = obj)
         IF ( .NOT. ASSOCIATED(obj) )     THEN
            self => NULL() 
         END IF      
      END SUBROUTINE releaseFTLinkedList
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION FTLinkedListDescription(self)  
         IMPLICIT NONE  
         CLASS(FTLinkedList)                         :: self
         CLASS(FTLinkedListRecord), POINTER          :: listRecord => NULL()
         CHARACTER(LEN=DESCRIPTION_CHARACTER_LENGTH) :: FTLinkedListDescription
         
         
         FTLinkedListDescription = ""
         IF(.NOT.ASSOCIATED(self % head)) RETURN
         
         listRecord              => self % head
         FTLinkedListDescription = TRIM(listRecord % recordObject % description())
         listRecord              => listRecord % next

         DO WHILE (ASSOCIATED(listRecord))
            FTLinkedListDescription = TRIM(FTLinkedListDescription) // &
                                       CHAR(13) // &
                                       TRIM(listRecord % recordObject % description())
            listRecord => listRecord % next
         END DO
      END FUNCTION FTLinkedListDescription    
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE printFTLinkedListDescription(self,iUnit)  
         IMPLICIT NONE  
         CLASS(FTLinkedList)                 :: self
         INTEGER                             :: iUnit
         CLASS(FTLinkedListRecord), POINTER  :: listRecord => NULL()
         LOGICAL                             :: circular
         
         
         IF(.NOT.ASSOCIATED(self % head)) RETURN
         
         IF(self % isCircular_) circular = .TRUE.
         CALL self % makeCircular(.FALSE.)
         
         listRecord => self % head

         DO WHILE (ASSOCIATED(listRecord))
            CALL listRecord % printDescription(iUnit)
            listRecord => listRecord % next
         END DO
         
         IF(circular) CALL self % makeCircular (.TRUE.)
         
      END SUBROUTINE printFTLinkedListDescription
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE reverseLinkedList(self)
!
!     ------------------------
!     Reverses the linked list
!     ------------------------
!
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(FTLinkedList) :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         CLASS(FTLinkedListRecord), POINTER :: current => NULL(), tmp => NULL()
         
         IF(.NOT.ASSOCIATED(self % head)) RETURN
         
         IF ( self % isCircular_ )     THEN
            self % head % previous => NULL()
            self % tail % next     => NULL() 
         END IF
         
         current  => self % head

         DO WHILE (ASSOCIATED(current))
            tmp                => current % next
            current % next     => current % previous
            current % previous => tmp
            current            => tmp
         END DO
         
         tmp => self % head
         self % head => self % tail
         self % tail => tmp
         
         IF ( self % isCircular_ )     THEN
            CALL self % makeCircular(.TRUE.) 
         END IF 
         
      END SUBROUTINE reverseLinkedList
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION allLinkedListObjects(self)  RESULT(array)
         IMPLICIT NONE  
!
!        ---------
!        Arguments
!        ---------
!
         CLASS (FTLinkedList)                 :: self
         CLASS(FTMutableObjectArray), POINTER :: array
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                            :: N
         CLASS(FTLinkedListRecord), POINTER :: listRecord => NULL()
         CLASS(FTObject)          , POINTER :: obj        => NULL()
         
         array => NULL()
         N = self % count()
         IF(N==0)     RETURN
         
         ALLOCATE(array)
         CALL array % initWithSize(arraySize  = N)
         
         listRecord => self % head

         DO WHILE (ASSOCIATED(listRecord))
            obj => listRecord % recordObject
            CALL array % addObject(obj)
            listRecord => listRecord % next
         END DO
         
      END FUNCTION allLinkedListObjects
!
!//////////////////////////////////////////////////////////////////////// 
! 
!      -----------------------------------------------------------------
!> Class name returns a string with the name of the type of the object
!>
!>  ### Usage:
!>
!>        PRINT *,  obj % className()
!>        if( obj % className = "FTLinkedList")
!>
      FUNCTION linkedListClassName(self)  RESULT(s)
         IMPLICIT NONE  
         CLASS(FTLinkedList)                        :: self
         CHARACTER(LEN=CLASS_NAME_CHARACTER_LENGTH) :: s
         
         s = "FTLinkedList"
 
      END FUNCTION linkedListClassName
!@mark -
! type conversions
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE castObjectToLinkedList(obj,cast) 
!
!     -----------------------------------------------------
!     Cast the base class FTObject to the LinkedList class
!     -----------------------------------------------------
!
         IMPLICIT NONE  
         CLASS(FTObject)    , POINTER :: obj
         CLASS(FTLinkedList), POINTER :: cast
         
         cast => NULL()
         SELECT TYPE (e => obj)
            TYPE is (FTLinkedList)
               cast => e
            CLASS DEFAULT
               
         END SELECT
         
      END SUBROUTINE castObjectToLinkedList
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION linkedListFromObject(obj) RESULT(cast)
!
!     -----------------------------------------------------
!     Cast the base class FTObject to the LinkedList class
!     -----------------------------------------------------
!
         IMPLICIT NONE  
         CLASS(FTObject)    , POINTER :: obj 
         CLASS(FTLinkedList), POINTER :: cast
         
         cast => NULL()
         SELECT TYPE (e => obj)
            TYPE is (FTLinkedList)
               cast => e
            CLASS DEFAULT
               
         END SELECT
         
      END FUNCTION linkedListFromObject
!
      END MODULE FTLinkedListClass
!
!@mark -
!
!//////////////////////////////////////////////////////////////////////// 
! 
!>An object for stepping through a linked list.
!>
!>###Definition (Subclass of FTObject):
!>   TYPE(FTLinkedListIterator) :: list
!>
!>
!>###Initialization
!>
!>         CLASS(FTLinkedList)        , POINTER :: list
!>         CLASS(FTLinkedListIterator), POINTER :: iterator
!>         ALLOCATE(iterator)
!>         CALL iterator % initWithFTLinkedList(list)
!>
!>###Accessors
!>
!>         ptr => iterator % list()
!>         ptr => iterator % object()
!>         ptr => iterator % currentRecord()
!>
!>###Iterating
!>
!>         CLASS(FTObject), POINTER :: objectPtr
!>         CALL iterator % setToStart
!>         DO WHILE (.NOT.iterator % isAtEnd())
!>            objectPtr => iterator % object()        ! if the object is wanted
!>            recordPtr => iterator % currentRecord() ! if the record is wanted
!>            
!>             Do something with object or record
!>
!>            CALL iterator % moveToNext() ! DON'T FORGET THIS!!
!>         END DO
!>
!>###Destruction
!>   
!>         CALL iterator % destruct() [Non Pointers]
!>         CALL release(iterator) [Pointers]
!
!//////////////////////////////////////////////////////////////////////// 
! 
      Module FTLinkedListIteratorClass
      USE FTLinkedListClass
      IMPLICIT NONE
!
!     -----------------
!     Class object type
!     -----------------
!
      TYPE, EXTENDS(FTObject) :: FTLinkedListIterator
         CLASS(FTLinkedList)      , POINTER :: list    => NULL()
         CLASS(FTLinkedListRecord), POINTER :: current => NULL()
!
!        ========         
         CONTAINS
!        ========
!
         PROCEDURE :: init           => initEmpty
         PROCEDURE :: initWithFTLinkedList
         PROCEDURE :: destruct       => destructIterator
         PROCEDURE :: isAtEnd        => FTLinkedListIsAtEnd
         PROCEDURE :: object         => FTLinkedListObject
         PROCEDURE :: currentRecord  => FTLinkedListCurrentRecord
         PROCEDURE :: linkedList     => returnLinkedList
         PROCEDURE :: className      => linkedListIteratorClassName
         PROCEDURE :: setLinkedList
         PROCEDURE :: setToStart
         PROCEDURE :: moveToNext
         PROCEDURE :: removeCurrentRecord
      END TYPE FTLinkedListIterator
      
      INTERFACE release
         MODULE PROCEDURE releaseFTLinkedListIterator 
      END INTERFACE  
!
!     ----------
!     Procedures
!     ----------
!
      CONTAINS 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE initEmpty(self) 
         IMPLICIT NONE 
         CLASS(FTLinkedListIterator)  :: self
!
!        --------------------------------------------
!        Always call the superclass initializer first
!        --------------------------------------------
!
         CALL self % FTObject % init()
!
!        ----------------------------------------------
!        Then call the initializations for the subclass
!        ----------------------------------------------
!
         self % list    => NULL()
         self % current => NULL()
         
      END SUBROUTINE initEmpty   
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE initWithFTLinkedList(self,list) 
         IMPLICIT NONE 
         CLASS(FTLinkedListIterator)  :: self
         CLASS(FTLinkedList), POINTER :: list
!
!        --------------------------------------------
!        Always call the superclass initializer first
!        --------------------------------------------
!
         CALL self % FTObject % init()
!
!        ----------------------------------------------
!        Then call the initializations for the subclass
!        ----------------------------------------------
!
         self % list    => NULL()
         self % current => NULL()
         CALL self % setLinkedList(list)
         CALL self % setToStart()
         
      END SUBROUTINE initWithFTLinkedList   
!
!////////////////////////////////////////////////////////////////////////
!
!< The destructor must not be called except at the end of destructors of
! subclasses.
!
      SUBROUTINE destructIterator(self)
          IMPLICIT NONE 
          CLASS(FTLinkedListIterator) :: self
          CLASS(FTObject), POINTER    :: obj
          
          CALL releaseMemberList(self)
          self % current => NULL()
!
!        ------------------------------------------
!        Always call the superclass destructor here
!        at the end of the subclass destructor.
!        ------------------------------------------
!
          CALL self % FTObject % destruct()
          
      END SUBROUTINE destructIterator
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE releaseMemberList(self)  
          IMPLICIT NONE  
          CLASS(FTLinkedListIterator) :: self
          CLASS(FTObject), POINTER    :: obj
          
          IF ( ASSOCIATED(self % list) )     THEN
             obj => self % list
             CALL releaseFTObject(self = obj)
             IF(.NOT. ASSOCIATED(obj)) self % list => NULL()
          END IF 
      END SUBROUTINE releaseMemberList
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
      SUBROUTINE releaseFTLinkedListIterator(self)  
         IMPLICIT NONE
         CLASS(FTLinkedListIterator) , POINTER :: self
         CLASS(FTObject)             , POINTER :: obj
         
         IF(.NOT. ASSOCIATED(self)) RETURN
         
         obj => self
         CALL releaseFTObject(self = obj)
         IF ( .NOT. ASSOCIATED(obj) )     THEN
            self => NULL() 
         END IF      
      END SUBROUTINE releaseFTLinkedListIterator
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE setToStart(self) 
         IMPLICIT NONE 
         CLASS(FTLinkedListIterator)  :: self
         self % current => self % list % head
      END SUBROUTINE setToStart 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE moveToNext(self) 
         IMPLICIT NONE 
         CLASS(FTLinkedListIterator)  :: self
         
         IF ( ASSOCIATED(self % current) )     THEN
            self % current => self % current % next
         ELSE 
            self % current => NULL() 
         END IF 
         
         IF ( ASSOCIATED(self % current, self % list % head) )     THEN
            self % current => NULL() 
         END IF 
         
      END SUBROUTINE moveToNext 
!
!////////////////////////////////////////////////////////////////////////
!
      LOGICAL FUNCTION FTLinkedListIsAtEnd(self)
         IMPLICIT NONE 
         CLASS(FTLinkedListIterator)  :: self
         IF ( ASSOCIATED(self % current) )     THEN
            FTLinkedListIsAtEnd = .false.
         ELSE
            FTLinkedListIsAtEnd = .true.
         END IF
      END FUNCTION FTLinkedListIsAtEnd   
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE setLinkedList(self,list)
         IMPLICIT NONE 
         CLASS(FTLinkedListIterator)  :: self
         CLASS(FTLinkedList), POINTER :: list
!
!        -----------------------------------
!        Remove current list if there is one
!        -----------------------------------
!
         IF ( ASSOCIATED(list) )     THEN
         
            IF ( ASSOCIATED(self % list, list) )     THEN
               CALL self % setToStart()
            ELSE IF( ASSOCIATED(self % list) )     THEN
               CALL releaseMemberList(self)
               self % list => list
               CALL self % list % retain()
               CALL self % setToStart
            ELSE
               self % list => list
               CALL self % list % retain()
               CALL self % setToStart()
            END IF 
            
         ELSE
         
            IF( ASSOCIATED(self % list) )     THEN
               CALL releaseMemberList(self)
            END IF 
            self % list => NULL()
            
         END IF
         
      END SUBROUTINE setLinkedList   
!
!////////////////////////////////////////////////////////////////////////
!
      FUNCTION returnLinkedList(self) RESULT(o)
         IMPLICIT NONE 
         CLASS(FTLinkedListIterator)  :: self
         CLASS(FTLinkedList), POINTER     :: o
         o => self % list
      END FUNCTION returnLinkedList 
!
!////////////////////////////////////////////////////////////////////////
!
      FUNCTION FTLinkedListObject(self) RESULT(o)
         IMPLICIT NONE 
         CLASS(FTLinkedListIterator)  :: self
         CLASS(FTObject), POINTER     :: o
         o => self % current % recordObject
      END FUNCTION FTLinkedListObject 
!
!////////////////////////////////////////////////////////////////////////
!
      FUNCTION FTLinkedListCurrentRecord(self) RESULT(o)
         IMPLICIT NONE 
         CLASS(FTLinkedListIterator)        :: self
         CLASS(FTLinkedListRecord), POINTER :: o
         o => self % current
      END FUNCTION FTLinkedListCurrentRecord 
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE removeCurrentRecord(self)  
         IMPLICIT NONE  
         CLASS(FTLinkedListIterator)        :: self
         CLASS(FTLinkedListRecord), POINTER :: r, n
         r => self % current
         n => self % current % next
         
         CALL self % list % removeRecord(r)
         self % current => n
         
      END SUBROUTINE removeCurrentRecord
!
!//////////////////////////////////////////////////////////////////////// 
! 
!      -----------------------------------------------------------------
!> Class name returns a string with the name of the type of the object
!>
!>  ### Usage:
!>
!>        PRINT *,  obj % className()
!>        if( obj % className = "FTLinkedListIterator")
!>
      FUNCTION linkedListIteratorClassName(self)  RESULT(s)
         IMPLICIT NONE  
         CLASS(FTLinkedListIterator)                :: self
         CHARACTER(LEN=CLASS_NAME_CHARACTER_LENGTH) :: s
         
         s = "FTLinkedListIterator"
 
      END FUNCTION linkedListIteratorClassName
      
      END MODULE FTLinkedListIteratorClass   