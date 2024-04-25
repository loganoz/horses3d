!
!////////////////////////////////////////////////////////////////////////
!
!      FTExceptionClass.f90
!      Created: January 29, 2013 5:06 PM 
!      By: David Kopriva  
!
!
!>An FTException object gives a way to pass generic
!>information about an exceptional situation.
!>
!>An FTException object gives a way to pass generic
!>information about an exceptional situation. Methods for
!>dealing with exceptions are defined in the SharedExceptionManagerModule
!>module.
!>
!>An FTException object wraps:
!>
!>- A severity indicator
!>- A name for the exception
!>- An optional dictionary that contains whatever information is deemed necessary.
!>
!>It is expected that classes will define exceptions that use instances
!>of the FTException Class.
!>
!>### Defined constants:
!>
!>-   FT_ERROR_NONE    = 0
!>-   FT_ERROR_WARNING = 1
!>-   FT_ERROR_FATAL   = 2
!>
!>### Initialization
!>
!>            CALL e  %  initFTException(severity,exceptionName,infoDictionary)
!>
!>Plus the convenience initializers, which automatically create a FTValueDictionary with a single key called "message":
!>
!>        CALL e % initWarningException(msg = "message")
!>        CALL e % initFatalException(msg = "message")
!>
!>Plus an assertion exception
!>
!>        CALL e % initAssertionFailureException(msg,expectedValueObject,observedValueObject,level)
!>
!>### Destruction
!>
!>        CALL e  %  release()
!>
!>###Setting the infoDictionary
!>
!>        CALL e  %  setInfoDictionary(infoDictionary)
!>###Getting the infoDictionary
!>
!>        dict => e % infoDictionary
!>###Getting the name of the exception
!>
!>        name = e % exceptionName()
!>###Getting the severity level of the exception
!>
!>        level = e % severity()
!> Severity levels are FT_ERROR_WARNING or FT_ERROR_FATAL
!>###Printing the exception
!>
!>        CALL e % printDescription()
!>
!
!////////////////////////////////////////////////////////////////////////
!
      Module FTExceptionClass
      USE FTStackClass
      USE FTDictionaryClass
      USE FTValueDictionaryClass
      USE FTLinkedListIteratorClass
      IMPLICIT NONE
!
!     ----------------
!     Global constants
!     ----------------
!
      INTEGER, PARAMETER :: FT_ERROR_NONE = 0, FT_ERROR_WARNING = 1, FT_ERROR_FATAL = 2
      INTEGER, PARAMETER :: ERROR_MSG_STRING_LENGTH = 132
      
      CHARACTER(LEN=21), PARAMETER :: FTFatalErrorException       = "FTFatalErrorException"
      CHARACTER(LEN=23), PARAMETER :: FTWarningErrorException     = "FTWarningErrorException"
      CHARACTER(LEN=27), PARAMETER :: FTAssertionFailureException = "FTAssertionFailureException"
!
!     ---------------
!     Error base type
!     ---------------
!
      TYPE, EXTENDS(FTObject) :: FTException
         INTEGER, PRIVATE                                :: severity_
         CHARACTER(LEN=ERROR_MSG_STRING_LENGTH), PRIVATE :: exceptionName_
         CLASS(FTDictionary), POINTER, PRIVATE           :: infoDictionary_ => NULL()
!
!        --------         
         CONTAINS
!        --------         
!
         PROCEDURE :: initFTException
         PROCEDURE :: initWarningException
         PROCEDURE :: initFatalException
         PROCEDURE :: initAssertionFailureException
         PROCEDURE :: destruct => destructException
         PROCEDURE :: setInfoDictionary
         PROCEDURE :: infoDictionary
         PROCEDURE :: exceptionName
         PROCEDURE :: severity
         PROCEDURE :: printDescription => printFTExceptionDescription
         PROCEDURE :: className => exceptionClassName
      END TYPE FTException
            
      INTERFACE cast
         MODULE PROCEDURE castToException
      END INTERFACE cast
      
      INTERFACE release
         MODULE PROCEDURE releaseFTException 
      END INTERFACE  
!
!     ========      
      CONTAINS
!     ========
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE initWarningException(self,msg)  
!
! ---------------------------------------------
!>A convenience initializer for a warning error 
!>that includes the key "message" in the
!>infoDictionary. Use this initializer as an 
!>example of how to write one's own exception.
! --------------------------------------------
!
         IMPLICIT NONE
         CLASS(FTException)                     :: self
         CHARACTER(LEN=*)                       :: msg
         
         CLASS(FTValueDictionary), POINTER :: userDictionary => NULL()
         CLASS(FTDictionary)     , POINTER :: dictPtr        => NULL()
            
         ALLOCATE(userDictionary)
         CALL userDictionary % initWithSize(64)
         CALL userDictionary % addValueForKey(msg,"message")
         
         dictPtr => userDictionary
         CALL self % initFTException(severity       = FT_ERROR_WARNING,&
                                     exceptionName  = FTWarningErrorException,&
                                     infoDictionary = dictPtr)
         CALL releaseMemberDictionary(self)
         
      END SUBROUTINE initWarningException
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE initFatalException(self,msg)  
!
! ---------------------------------------------
!>A convenience initializer for a fatal error 
!>that includes the key "message" in the
!>infoDictionary.Use this initializer as an 
!>example of how to write one's own exception.
! --------------------------------------------
!
         IMPLICIT NONE
         CLASS(FTException)                     :: self
         CHARACTER(LEN=*)                       :: msg
         
         CLASS(FTValueDictionary), POINTER :: userDictionary => NULL()
         CLASS(FTDictionary)     , POINTER :: dictPtr        => NULL()
            
         ALLOCATE(userDictionary)
         CALL userDictionary % initWithSize(8)
         CALL userDictionary % addValueForKey(msg,"message")
         
         dictPtr => userDictionary
         CALL self % initFTException(severity       = FT_ERROR_FATAL,&
                                     exceptionName  = FTFatalErrorException,&
                                     infoDictionary = dictPtr)
         
         CALL releaseMemberDictionary(self)
         
      END SUBROUTINE initFatalException
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE initFTException(self,severity,exceptionName,infoDictionary)
!
! -----------------------------------
!>The main initializer for the class 
! -----------------------------------
!
         IMPLICIT NONE
         CLASS(FTException)                     :: self
         INTEGER                                :: severity
         CHARACTER(LEN=*)                       :: exceptionName
         CLASS(FTDictionary), POINTER, OPTIONAL :: infoDictionary
         
         CALL self  %  FTObject  %  init()
         
         self  %  severity_        = severity
         self  %  exceptionName_   = exceptionName
         self  %  infoDictionary_  => NULL()
         IF(PRESENT(infoDictionary) .AND. ASSOCIATED(infoDictionary))   THEN 
            CALL self % setInfoDictionary(infoDictionary)
         END IF 
         
      END SUBROUTINE initFTException
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE initAssertionFailureException(self,msg,expectedValueObject,observedValueObject,level)
!
! ------------------------------------------------
!>A convenience initializer for an assertion error 
!>that includes the keys:
!>
!>-"message"
!>-"expectedValue"
!>-"observedValue"
!>
!>in the infoDictionary
!
! ------------------------------------------------
!
         IMPLICIT NONE
         CLASS(FTException)      :: self
         CLASS(FTValue), POINTER :: expectedValueObject, ObservedValueObject
         INTEGER                 :: level
         CHARACTER(LEN=*)        :: msg
         
         CLASS(FTValueDictionary), POINTER :: userDictionary => NULL()
         CLASS(FTDictionary)     , POINTER :: dictPtr        => NULL()
         CLASS(FTObject)         , POINTER :: objectPtr      => NULL()
            
         ALLOCATE(userDictionary)
         CALL userDictionary % initWithSize(8)
         CALL userDictionary % addValueForKey(msg,"message")
         objectPtr => expectedValueObject
         CALL userDictionary % addObjectForKey(object = objectPtr,key = "expectedValue")
         objectPtr => ObservedValueObject
         CALL userDictionary % addObjectForKey(object = objectPtr,key = "observedValue")
         
         dictPtr => userDictionary
         CALL self % initFTException(severity       = level,&
                                     exceptionName  = FTAssertionFailureException,&
                                     infoDictionary = dictPtr)
         
         CALL releaseMemberDictionary(self)
         
      END SUBROUTINE initAssertionFailureException
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE destructException(self)
!
! -------------------------------------------------------------
!>The destructor for the class. Do not call this directly. Call
!>the release() procedure instead
! -------------------------------------------------------------
!

         IMPLICIT NONE  
         CLASS(FTException)       :: self

         CALL releaseMemberDictionary(self)
         CALL self % FTObject % destruct()
         
      END SUBROUTINE destructException 
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
      SUBROUTINE releaseFTException(self)  
         IMPLICIT NONE
         CLASS(FTException) , POINTER :: self
         CLASS(FTObject)    , POINTER :: obj
         
         IF(.NOT. ASSOCIATED(self)) RETURN
         
         obj => self
         CALL releaseFTObject(self = obj)
         IF ( .NOT. ASSOCIATED(obj) )     THEN
            self => NULL() 
         END IF      
      END SUBROUTINE releaseFTException
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE setInfoDictionary( self, dict )  
!
! ---------------------------------------------
!>Sets and retains the exception infoDictionary
! ---------------------------------------------
!
         IMPLICIT NONE
         CLASS(FTException)           :: self
         CLASS(FTDictionary), POINTER :: dict
         
         IF(ASSOCIATED(self % infoDictionary_)) CALL releaseMemberDictionary(self)
         self  %  infoDictionary_ => dict
         CALL self  %  infoDictionary_  %  retain()
      END SUBROUTINE setInfoDictionary
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE releaseMemberDictionary(self)  
         IMPLICIT NONE  
         CLASS(FTException)       :: self
         CLASS(FTObject), POINTER :: obj
         
         IF(ASSOCIATED(self % infoDictionary_))   THEN
            obj => self % infoDictionary_
            CALL releaseFTObject(self = obj)
            IF(.NOT. ASSOCIATED(obj)) self% infoDictionary_ => NULL()
         END IF
      END SUBROUTINE releaseMemberDictionary
!
!//////////////////////////////////////////////////////////////////////// 
! 
     FUNCTION infoDictionary(self)
!
! ---------------------------------------------
!>Returns the exception's infoDictionary. Does
!>not transfer ownership/reference count is 
!>unchanged.
! ---------------------------------------------
!
        IMPLICIT NONE  
        CLASS(FTException) :: self
        CLASS(FTDictionary), POINTER :: infoDictionary
        
        infoDictionary => self % infoDictionary_
        
     END FUNCTION infoDictionary
!
!//////////////////////////////////////////////////////////////////////// 
! 
     FUNCTION exceptionName(self)  
!
! ---------------------------------------------
!>Returns the string representing the name set
!>for the exception.
! ---------------------------------------------
!
        IMPLICIT NONE  
        CLASS(FTException) :: self
        CHARACTER(LEN=ERROR_MSG_STRING_LENGTH) :: exceptionName
        exceptionName = self % exceptionName_
     END FUNCTION exceptionName
!
!//////////////////////////////////////////////////////////////////////// 
! 
     INTEGER FUNCTION severity(self)  
!
! ---------------------------------------------
!>Returns the severity level of the exception.
! ---------------------------------------------
!
        IMPLICIT NONE  
        CLASS(FTException) :: self
        severity = self % severity_
     END FUNCTION severity    
!
!//////////////////////////////////////////////////////////////////////// 
! 
     SUBROUTINE printFTExceptionDescription(self,iUnit)  
!
! ----------------------------------------------
!>A basic printing of the exception and the info
!>held in the infoDicitonary.
! ----------------------------------------------
!
        IMPLICIT NONE  
        CLASS(FTException) :: self
        INTEGER            :: iUnit
        
        CLASS(FTDictionary), POINTER :: dict => NULL()
        
!        WRITE(iUnit,*) "-------------"
        WRITE(iUnit,*) " "
        WRITE(iUnit,*) "Exception Named: ", TRIM(self  %  exceptionName())
        dict => self % infoDictionary()
        IF(ASSOCIATED(dict)) CALL dict % printDescription(iUnit)
        
     END SUBROUTINE printFTExceptionDescription     
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE castToException(obj,cast) 
!
! -----------------------------------------------------
!>Cast the base class FTObject to the FTException class
! -----------------------------------------------------
!
         IMPLICIT NONE  
         CLASS(FTObject)   , POINTER :: obj
         CLASS(FTException), POINTER :: cast
         
         cast => NULL()
         SELECT TYPE (e => obj)
            TYPE is (FTException)
               cast => e
            CLASS DEFAULT
               
         END SELECT
         
      END SUBROUTINE castToException
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION exceptionFromObject(obj) RESULT(cast)
!
!     -----------------------------------------------------
!     Cast the base class FTObject to the FTException class
!     -----------------------------------------------------
!
         IMPLICIT NONE  
         CLASS(FTObject)   , POINTER :: obj
         CLASS(FTException), POINTER :: cast
         
         cast => NULL()
         SELECT TYPE (e => obj)
            TYPE is (FTException)
               cast => e
            CLASS DEFAULT
               
         END SELECT
         
      END FUNCTION exceptionFromObject
!
!//////////////////////////////////////////////////////////////////////// 
! 
!      -----------------------------------------------------------------
!> Class name returns a string with the name of the type of the object
!>
!>  ### Usage:
!>
!>        PRINT *,  obj % className()
!>        if( obj % className = "FTException")
!>
      FUNCTION exceptionClassName(self)  RESULT(s)
         IMPLICIT NONE  
         CLASS(FTException)                         :: self
         CHARACTER(LEN=CLASS_NAME_CHARACTER_LENGTH) :: s
         
         s = "FTException"
 
      END FUNCTION exceptionClassName

      END Module FTExceptionClass
!
!//////////////////////////////////////////////////////////////////////// 
! 
!@mark -
     
      Module SharedExceptionManagerModule
!>
!>All exceptions are posted to the SharedExceptionManagerModule. 
!>
!>To use exceptions,first initialize it
!>        CALL initializeFTExceptions
!>From that point on, all exceptions will be posted there. Note that the
!>FTTestSuiteManager class will initialize the SharedExceptionManagerModule,
!>so there is no need to do the initialization separately if the FTTestSuiteManager
!>class has been initialized.
!>
!>The exceptions are posted to a stack. To access the exceptions they will be
!>peeked or popped from that stack.
!>
!>###Initialization
!>        CALL initializeFTExceptions
!>###Finalization
!>        CALL destructFTExceptions
!>###Throwing an exception
!>         CALL throw(exception)
!>###Getting the number of exceptions
!>         n = errorCount()
!>###Catching all exceptions
!>         IF(catch())     THEN
!>            Do something with the exceptions
!>         END IF
!>###Getting the named exception caught
!>         CLASS(FTException), POINTER :: e
!>         e => errorObject()
!>###Popping the top exception
!>         e => popLastException()
!>###Peeking the top exception
!>         e => peekLastException()
!>###Catching an exception with a given name
!>         IF(catch(name))   THEN
!>            !Do something with the exception, e.g.
!>            e              => errorObject()
!>            d              => e % infoDictionary()
!>            userDictionary => valueDictionaryFromDictionary(dict = d)
!>            msg = userDictionary % stringValueForKey("message",FTDICT_KWD_STRING_LENGTH)
!>         END IF
!>###Printing all exceptions
!>      call printAllExceptions
!>         
      USE FTExceptionClass
      IMPLICIT NONE  
!
!     --------------------
!     Global error stack  
!     --------------------
!
      TYPE(FTStack)    , POINTER, PRIVATE :: errorStack    => NULL()
      TYPE(FTException), POINTER, PRIVATE :: currentError_ => NULL()
      
      INTERFACE catch
         MODULE PROCEDURE catchAll
         MODULE PROCEDURE catchErrorWithName
      END INTERFACE catch
      
      PRIVATE :: catchAll, catchErrorWithName
!
!     ========      
      CONTAINS
!     ========
!
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE initializeFTExceptions
!
!>Called at start of execution. Will be called automatically if an 
!>exception is thrown.
!
         IMPLICIT NONE
         IF ( .NOT.ASSOCIATED(errorStack) )     THEN
            ALLOCATE(errorStack)
            CALL errorStack % init()
            currentError_ => NULL()
         END IF 
      END SUBROUTINE initializeFTExceptions
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE destructFTExceptions
!
!>Called at the end of execution. This procedure will announce if there
!>are uncaught exceptions raised and print them.
!
         IMPLICIT NONE
         CLASS(FTObject), POINTER :: obj
!  
!        --------------------------------------------------
!        First see if there are any uncaught exceptions and
!        report them if there are.
!        --------------------------------------------------
!
         IF ( catch() )     THEN
           PRINT *
           PRINT *,"   ***********************************"
           IF(errorStack % COUNT() == 1)     THEN
              PRINT *, "   An uncaught exception was raised:"
           ELSE
              PRINT *, "   Uncaught exceptions were raised:"
           END IF
           PRINT *,"   ***********************************"
           PRINT *
           !DEBUG CALL printAllExceptions
           CALL errorStack % printDescription(iUnit = 6)!DEBUG
         END IF 
!
!        -----------------------
!        Destruct the exceptions
!        -----------------------
!

          obj => errorStack
          CALL releaseFTObject(self = obj)
          IF(.NOT. ASSOCIATED(obj)) errorStack => NULL()
          CALL releaseCurrentError
        
      END SUBROUTINE destructFTExceptions
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE throw(exceptionToThrow)
!
!>Throws the exception: exceptionToThrow
!
         IMPLICIT NONE  
         CLASS(FTException), POINTER :: exceptionToThrow
         CLASS(FTObject)   , POINTER :: ptr => NULL()
         
         IF ( .NOT.ASSOCIATED(errorStack) )     THEN
            CALL initializeFTExceptions 
         END IF 
         
         ptr => exceptionToThrow
         CALL errorStack % push(ptr)
         
      END SUBROUTINE throw
!
!//////////////////////////////////////////////////////////////////////// 
! 
      LOGICAL FUNCTION catchAll()
!
! -------------------------------------------
!>Returns .TRUE. if there are any exceptions.
! -------------------------------------------
!
         IMPLICIT NONE
         
         IF ( .NOT.ASSOCIATED(errorStack) )     THEN
            catchAll = .FALSE.
            RETURN 
         END IF 
         
         catchAll = .false.
         IF ( errorStack % count() > 0 )     THEN
            catchAll = .true.
         END IF
         CALL releaseCurrentError
         
      END FUNCTION catchAll
!
!//////////////////////////////////////////////////////////////////////// 
! 
      INTEGER FUNCTION errorCount()
!
! ------------------------------------------
!>Returns the number of exceptions that have 
!>been thrown.
! ------------------------------------------
!
         IMPLICIT NONE
                  
         IF ( .NOT.ASSOCIATED(errorStack) )     THEN
            CALL initializeFTExceptions 
         END IF 

         errorCount = errorStack % count() 
      END FUNCTION    
!
!//////////////////////////////////////////////////////////////////////// 
! 
      LOGICAL FUNCTION catchErrorWithName(exceptionName)
!
! --------------------------------------------
!>Returns .TRUE. if there is an exception with
!>the requested name. If so, it pops the 
!>exception and saves the pointer to it so that
!>it can be accessed with the currentError()
!>function.
! --------------------------------------------
!
     
         IMPLICIT NONE  
         CHARACTER(LEN=*) :: exceptionName
         
         TYPE(FTLinkedListIterator)   :: iterator
         CLASS(FTLinkedList), POINTER :: ptr => NULL()
         CLASS(FTObject)    , POINTER :: obj => NULL()
         CLASS(FTException) , POINTER :: e   => NULL()
         
         catchErrorWithName = .false.
                  
         IF ( .NOT.ASSOCIATED(errorStack) )     THEN
            CALL initializeFTExceptions 
            RETURN 
         END IF 
         
         IF ( errorStack % COUNT() == 0 )     THEN
            RETURN 
         END IF 

         ptr => errorStack
         CALL iterator % initWithFTLinkedList(ptr)
         CALL iterator % setToStart()
         
         DO WHILE (.NOT.iterator % isAtEnd())
            obj => iterator % object()
            CALL cast(obj,e)
            IF ( e % exceptionName() == exceptionName )     THEN
               CALL setCurrentError(e)
               catchErrorWithName = .true.
               CALL errorStack % remove(obj)
               EXIT
           END IF 
           CALL iterator % moveToNext()
         END DO
         
         CALL iterator % destruct()
         
      END FUNCTION catchErrorWithName
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION errorObject()
!
! -------------------------------------------
!>Returns a pointer to the current exception.
! -------------------------------------------
!
         IMPLICIT NONE
         CLASS(FTException), POINTER :: errorObject
         
         IF ( .NOT.ASSOCIATED(errorStack) )     THEN
            CALL initializeFTExceptions 
         END IF 
         
         errorObject => currentError_
      END FUNCTION errorObject
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE setCurrentError(e)  
         IMPLICIT NONE  
         CLASS(FTException) , POINTER :: e
!
!        --------------------------------------------------------------
!        Check first to see if there is a current error. Since it
!        is retained, the current one must be released before resetting
!        the pointer.
!        --------------------------------------------------------------
!
         CALL releaseCurrentError
!
!        ------------------------------------
!        Set the pointer and retain ownership
!        ------------------------------------
!
         currentError_ => e
         CALL currentError_ % retain()
         
      END SUBROUTINE setCurrentError
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION popLastException()
!
! ----------------------------------------------------------------
!>Get the last exception posted. This is popped from the stack.
!>The caller is responsible for releasing the object after popping
! ----------------------------------------------------------------
!
         IMPLICIT NONE  
         CLASS(FTException), POINTER :: popLastException
         CLASS(FTObject)   , POINTER :: obj => NULL()
         
         obj => NULL()
         popLastException => NULL()
         IF ( .NOT.ASSOCIATED(errorStack) )     THEN
            CALL initializeFTExceptions 
         ELSE
            CALL errorStack % pop(obj)
            IF(ASSOCIATED(obj)) CALL cast(obj,popLastException)
         END IF 
         
      END FUNCTION popLastException
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION peekLastException()  
!
! ----------------------------------------------------------------
!>Get the last exception posted. This is NOT popped from the stack.
!>The caller does not own the object.
! ----------------------------------------------------------------
!
         IMPLICIT NONE  
         CLASS(FTException), POINTER :: peekLastException
         CLASS(FTObject)   , POINTER :: obj => NULL()
         
         IF ( .NOT.ASSOCIATED(errorStack) )     THEN
            CALL initializeFTExceptions 
         END IF 
         
         peekLastException => NULL()
         obj => errorStack % peek()
         CALL cast(obj,peekLastException)
         
      END FUNCTION peekLastException
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE printAllExceptions  
         IMPLICIT NONE  
         TYPE(FTLinkedListIterator)   :: iterator
         CLASS(FTLinkedList), POINTER :: list      => NULL()
         CLASS(FTObject)    , POINTER :: objectPtr => NULL()
         CLASS(FTException) , POINTER :: e         => NULL()
           
        list => errorStack
        CALL iterator % initWithFTLinkedList(list)
!
!       ----------------------------------------------------
!       Write out the descriptions of each of the exceptions
!       ----------------------------------------------------
!
        CALL iterator % setToStart
        DO WHILE (.NOT.iterator % isAtEnd())
            objectPtr => iterator % object()
            CALL cast(objectPtr,e)
            CALL e % printDescription(6)
            CALL iterator % moveToNext()
         END DO
         
         CALL iterator % destruct() !iterator is not a pointer
            
      END SUBROUTINE printAllExceptions
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE releaseCurrentError
         IMPLICIT NONE
         CLASS(FTObject), POINTER :: obj
         
         IF ( ASSOCIATED(currentError_) )     THEN
           obj => currentError_
           CALL releaseFTObject(self = obj)
           IF(.NOT. ASSOCIATED(obj)) currentError_ => NULL()
         END IF 
 
      END SUBROUTINE releaseCurrentError

      END MODULE SharedExceptionManagerModule    
      