!
!////////////////////////////////////////////////////////////////////////
!
!
!
!>FTObject is the root class for all object types.
!>
!>Overview
!>--------
!>
!>FTObject defines the basic methods that are essential for reference counted objects.
!>
!>FTObject is generally not going to be instantiated by itself, but rather it will 
!>be subclassed and one will work with instances of the subclasses. 
!>Otherwise, pointers of type FTObject that point to instances of subclasses
!>will be stored in the container classes.
!>
!>
!>Tasks
!>-----
!>
!>
!>- init()
!>
!>     Initializes an object and any memory that it needs to allocate, etc. 
!>     Should be orrerrided in subclasses.The base class implementation does nothing but
!>     increase the reference count of the object.
!>
!>- destruct()
!>
!>     Destructor of the object, which releases and deallocates owned objects and memory.
!>     Should be overridden in subclasses. The base class implementation does nothing but
!>     decrease the reference count of the object.
!>
!>- printDescription(iUnit)
!>
!>     Prints a description of the object to a specified file unit. The base class implementation
!>     does nothing but print "FTObject"
!>
!>- copy()
!>
!>     Creates a copy (pointer) to the object of CLASS(FTObject) sourced with the object.
!>
!>- retain()
!>
!>     Increases the reference count of the object. Any procedure or object that retain()'s
!>     an object gains an ownership stake in that object. This procedure is not overridable.
!>
!>- release()
!>
!>     Decreases the reference count of an object. To be called only by objects or procedures
!>     that have ownership in an object pointer, i.e., for which init() or retain() have been called.
!>     Override this procedure in subclasses for releasing the actual type.
!>
!>- isUnreferenced()
!>
!>     Test to see if there are no more owners of an object.
!>
!>- refCount()
!>
!>     Returns the number of owners of an object. Usually this is of interest only for debugging purposes.
!>     This procedure is not overridable.
!>     
!>
!>Subclassing FTObject
!>--------------------
!>
!>In general, subclasses of FTObject override
!>
!>- init()
!>- destruct()
!>- printDescription()
!>- release()
!>
!>They should also provide a cast() subroutine to convert from the base class to a subclass.
!>The cast() routine can look something like
!>
!>     SUBROUTINE castToSubclass(obj,cast) 
!>        IMPLICIT NONE  
!>        CLASS(FTObject), POINTER :: obj
!>        CLASS(SubClass), POINTER :: cast
!>        
!>        cast => NULL()
!>        SELECT TYPE (e => obj)
!>           TYPE is (SubClass)
!>              cast => e
!>           CLASS DEFAULT
!>              
!>        END SELECT
!>        
!>     END SUBROUTINE castToSubclass
!>
!>
!>## Subclassing init
!>
!>The init() procedure performs subclass specific operations to initialize an object.
!>
!>Subclasses that override init() must include 
!>a call to the super class method. For example, overriding init looks like
!>
!>     SUBROUTINE initSubclass(self) 
!>        IMPLICIT NONE
!>        CLASS(Subclass) :: self
!>        
!>        CALL self % FTObject % init()
!>        Allocate and initialize all member objects
!>        ... Other Subclass specific code
!>     END SUBROUTINE initSubclass
!>
!>## Subclassing destruct
!>
!>The destruct() procedure reverses the operations done in the init() procedure. It releases and
!>deallocates any pointers that it owns.  Subclasses that override destruct() must include 
!>a call to the super class method. For example, overriding destruct looks like
!>
!>     SUBROUTINE destructSubclass(self) 
!>        IMPLICIT NONE
!>        CLASS(Subclass) :: self
!>        
!>        Release and deallocate (if necessary) all member objects
!>        CALL self % FTObject % destruct()
!>        
!>     END SUBROUTINE destructSubclass
!>
!>## Subclassing printDescription(iUnit)
!>
!>printDescription is a method whose existence is to support debugging. Call printDescription(iUnit)
!>on any objects owned by self for a cascading of what is stored in the object.
!>
!>
!>## Casting an object from the base to a subclass
!>
!>Container classes and the copy function return pointers to a CLASS(FTObject). To use
!>any subclass features one must "cast" to the subclass. We like to have a specific 
!>cast routine to do this as painlessly as possible. Each subclass should include a 
!>SUBROUTINE like this:
!>
!>     SUBROUTINE castToSubclass(obj,cast) 
!>        IMPLICIT NONE  
!>        CLASS(FTObject), POINTER :: obj
!>        CLASS(Subclass), POINTER :: cast
!>        cast => NULL()
!>        SELECT TYPE (e => obj)
!>           TYPE is (Subclass)
!>              cast => e
!>           CLASS DEFAULT
!>        END SELECT
!>     END SUBROUTINE castToValue
!>## Subclassing className
!>
!>The className() procedure returns the name of the class.
!>
!>Subclasses should override className() !>
!>
!>Created: January 7, 2013 11:30 AM 

!>@author 
!>David A. Kopriva 
!
      Module FTObjectClass 
      IMPLICIT NONE
!
!     ---------
!     Constants
!     ---------
!
      ! Until all compilers can do allocatables
      INTEGER, PARAMETER :: DESCRIPTION_CHARACTER_LENGTH = 1024
      INTEGER, PARAMETER :: CLASS_NAME_CHARACTER_LENGTH  = 32
!
!     --------------------------
!     Derived type for the class
!     --------------------------
!
      TYPE FTObject
         INTEGER, PRIVATE   :: refCount_
!
!        ========         
         CONTAINS
!        ========
!
         PROCEDURE                  :: init             => initFTObject
         PROCEDURE                  :: destruct         => destructFTObject
         PROCEDURE                  :: description      => FTObjectDescription
         PROCEDURE                  :: printDescription => printFTObjectDescription
         PROCEDURE                  :: className
         
         PROCEDURE, NON_OVERRIDABLE :: copy    => copyFTObject
         PROCEDURE, NON_OVERRIDABLE :: retain  => retainFTObject
         PROCEDURE, NON_OVERRIDABLE :: isUnreferenced
         PROCEDURE, NON_OVERRIDABLE :: refCount
      END TYPE FTObject
      
      PRIVATE :: copyFTObject
!
!     =====================
      CONTAINS ! Procedures
!     =====================
!
!
!////////////////////////////////////////////////////////////////////////
!
!
!     --------------------------------------------------------------
!>   Generic Name: init()
!>
!>   Initializes the object. The base class initialization does 
!>   nothing but set the reference count to one.
!     --------------------------------------------------------------
!
      SUBROUTINE initFTObject(self)
         IMPLICIT NONE 
         CLASS(FTObject) :: self
         self % refCount_   = 1
      END SUBROUTINE initFTObject
!
!////////////////////////////////////////////////////////////////////////
!
!
!     --------------------------------------------------------------------
!>   Generic Name: destruct()
!>
!>   The destructor for the class. The base class destructor does nothing.
!     --------------------------------------------------------------------
!
      SUBROUTINE destructFTObject(self)
         IMPLICIT NONE 
         CLASS(FTObject) :: self
    END SUBROUTINE destructFTObject
!
!////////////////////////////////////////////////////////////////////////
!
!
!      -----------------------------------------------------------------
!> Retain increases the reference count by one and implies ownership
!>  to the caller.
!>  ### Usage:
!>        CALL obj\ % retain()
!      -----------------------------------------------------------------
!
       SUBROUTINE retainFTObject(self)
         IMPLICIT NONE 
         CLASS(FTObject) :: self
         self % refCount_ = self % refCount_ + 1
      END SUBROUTINE retainFTObject
!
!////////////////////////////////////////////////////////////////////////
!
!
!      ---------------------------------------------------------------------------
!>     releaseFTObject decreases the reference count by one and implies 
!>     relinquishing ownership by the caller. Call this if control
!>     over the existence of an object pointer is no longer desired by the caller.
!>     When the reference count goes to zero, the destructor of the object
!>     is called automatically and the object is deallocated.
!      ---------------------------------------------------------------------------
!
       RECURSIVE SUBROUTINE releaseFTObject(self)
         IMPLICIT NONE 
         CLASS(FTObject), POINTER  :: self
         
         
         IF(.NOT. ASSOCIATED(self)) RETURN
         
         self % refCount_ = self % refCount_ - 1
         
         IF ( self % refCount_ < 0 )     THEN
            PRINT *, "Attempt to release object with refCount = 0"
            CALL self % printDescription(6)
            PRINT *, "--------------------------------------------"
            PRINT *, " "
            RETURN 
         END IF
         
         IF ( self % refCount_ == 0 )     THEN
            CALL self % destruct()
            DEALLOCATE(self)
            self => NULL()
         END IF 
          
      END SUBROUTINE releaseFTObject
!
!//////////////////////////////////////////////////////////////////////// 
! 
!      -----------------------------------------------------------------
!> Class name returns a string with the name of the type of the object
!>
!>  ### Usage:
!>
!>        PRINT *,  obj % className()
!>        if( obj % className = "FTObject")
!>
!      -----------------------------------------------------------------
!
      FUNCTION className(self)  RESULT(s)
         IMPLICIT NONE  
         CLASS(FTObject)                            :: self
         CHARACTER(LEN=CLASS_NAME_CHARACTER_LENGTH) :: s
         
         s = "FTObject"
 
      END FUNCTION className
!
!////////////////////////////////////////////////////////////////////////
!
!
!     -----------------------------------------------------------------------
!>    Owners of objects should call isUnreferenced after releasing a 
!>     pointer object. If true, the object should be deallocated and then
!>     set to point to NULL()
!>
!>     ### Usage: ###
!>
!>          IF ( v % isUnreferenced() )     THEN
!>             DEALLOCATE(v)
!>             v => NULL()
!>          END IF
!>    
!     -----------------------------------------------------------------------
!
      LOGICAL FUNCTION isUnreferenced(self)
         IMPLICIT NONE 
         CLASS(FTObject) :: self
         IF ( self % refCount_ == 0 )     THEN
            isUnreferenced = .true.
         ELSE
            isUnreferenced = .false.
         END IF
          
      END FUNCTION isUnreferenced
!
!////////////////////////////////////////////////////////////////////////
!
!
!     -----------------------------------------------------------------
!>   Returns the reference count for the object. Normally this is done
!>    only for debugging purposes.
!<
!     -----------------------------------------------------------------
!
      INTEGER FUNCTION refCount(self)
         IMPLICIT NONE 
         CLASS(FTObject) :: self
         refCount = self % refCount_
      END FUNCTION refCount 
!
!//////////////////////////////////////////////////////////////////////// 
! 
!
!     ----------------------------------------------------------------------
!>   Returns a character string of length DESCRIPTION_CHARACTER_LENGTH that
!>    represents the object. the base class implementation returns an empty
!>    string. Note that if the description is too long, the expected string
!>    will be truncated. In general, one wants to use printDescription.
!<
!     ----------------------------------------------------------------------
!
      FUNCTION FTObjectDescription(self)
         IMPLICIT NONE  
         CLASS(FTObject)    :: self
         CHARACTER(LEN=DESCRIPTION_CHARACTER_LENGTH) :: FTObjectDescription
         FTObjectDescription = " "
      END FUNCTION FTObjectDescription
!
!//////////////////////////////////////////////////////////////////////// 
! 
!
!     ------------------------------------------------------------------------------------
!>   Generic Name: printDescription()
!>
!>   Prints a string to unit iUnit that represents the contents of the object. FTObject's
!>    description simply prints its name. Override this in subclasses to print something
!>    useful. 
!<
!     ------------------------------------------------------------------------------------
!
      SUBROUTINE printFTObjectDescription(self,iUnit)
         IMPLICIT NONE  
         CLASS(FTObject) :: self
         INTEGER         :: iUnit
         WRITE(iUnit,*) "FTObject"
      END SUBROUTINE printFTObjectDescription
!
!//////////////////////////////////////////////////////////////////////// 
! 
!
!     --------------------------------------------------------------------
!>   Base class implementation of the assignment function. Call this from
!>    within any subclasses copy assignment function. All FTObject's 
!>    implementation does is set
!>    the reference count to one, implying no additional ownership to the 
!>    caller that is creating the copy.
!<
!     --------------------------------------------------------------------
!
      FUNCTION copyFTObject(self) RESULT(copy)
         IMPLICIT NONE  
         CLASS(FTObject), INTENT(IN) :: self
         CLASS(FTObject), POINTER    :: copy
         
         ALLOCATE(copy)
         CALL initFTObject(self = copy)
      END FUNCTION copyFTObject
      
      END MODULE FTObjectClass