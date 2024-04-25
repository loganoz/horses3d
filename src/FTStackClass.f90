!
!////////////////////////////////////////////////////////////////////////
!
!      FTStackClass.f90
!      Created: January 25, 2013 12:56 PM 
!      By: David Kopriva  
!
!>Inherits from FTLinkedListClass : FTObjectClass
!>
!>##Definition (Subclass of FTLinkedListClass):
!>   TYPE(FTStack) :: list
!>
!>#Usage:
!>
!>##Initialization
!>
!>      ALLOCATE(stack)  If stack is a pointer
!>      CALL stack  %  init()
!>
!>##Destruction
!>      CALL release(stack) [Pointers]
!>      CALL stack % destruct() [Non pointers]
!>
!>##Pushing an object onto the stack
!>
!>      TYPE(FTObject) :: objectPtr
!>      objectPtr => r1
!>      CALL stack % push(objectPtr)
!>
!>##Peeking at the top of the stack
!>
!>      objectPtr => stack % peek()  No change of ownership
!>      SELECT TYPE(objectPtr)
!>         TYPE is (*SubclassType*)
!>            … Do something with ObjectPtr as subclass
!>         CLASS DEFAULT
!>            … Problem with casting
!>      END SELECT
!>
!>##Popping the top of the stack
!>
!>      objectPtr => stack % pop()  Ownership transferred to caller
!>      SELECT TYPE(objectPtr)
!>         TYPE is (*SubclassType*)
!>            … Do something with ObjectPtr as subclass
!>         CLASS DEFAULT
!>            … Problem with casting
!>      END SELECT
!
!////////////////////////////////////////////////////////////////////////
!
      Module FTStackClass
      USE FTLinkedListClass
      IMPLICIT NONE
      
      TYPE, EXTENDS(FTLinkedList) :: FTStack
!
!        ========         
         CONTAINS
!        ========
!
         PROCEDURE :: init             => initFTStack
         PROCEDURE :: printDescription => printStackDescription
         PROCEDURE :: className        => stackClassName
         PROCEDURE :: push
         PROCEDURE :: pop
         PROCEDURE :: peek
      END TYPE FTStack
      
      INTERFACE release
         MODULE PROCEDURE  releaseFTStack
      END INTERFACE  
!
!     ----------
!     Procedures
!     ----------
!
!     ========
      CONTAINS
!     ========
!
!
!------------------------------------------------
!> Public, generic name: init()
!>
!> Initialize the stack.
!------------------------------------------------
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE initFTStack(self) 
         IMPLICIT NONE 
         CLASS(FTStack) :: self
!
!        --------------------------------------------
!        Call the initializer of the superclass first
!        --------------------------------------------
!
         CALL self % FTLinkedList % init()
!
!        ---------------------------------
!        Then initialize ivars of subclass
!        ---------------------------------
!
         !None to initialize
         
      END SUBROUTINE initFTStack
!
!------------------------------------------------------
!> Public, generic name: release(self)
!>
!> Call release(self) on an object to release control
!> of a pointer object. If its reference count is zero, 
!> then it is deallocated.
!------------------------------------------------------
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE releaseFTStack(self)  
         IMPLICIT NONE
         TYPE(FTStack)  , POINTER :: self
         CLASS(FTObject), POINTER :: obj
         
         IF(.NOT. ASSOCIATED(self)) RETURN
         
         obj => self
         CALL releaseFTObject(self = obj)
         IF ( .NOT. ASSOCIATED(obj) )     THEN
            self => NULL() 
         END IF      
      END SUBROUTINE releaseFTStack
!
!     -----------------------------------
!     push: Push an object onto the stack
!     -----------------------------------
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE push(self,obj)
!
!        ----------------------------------
!        Add object to the head of the list
!        ----------------------------------
!
         IMPLICIT NONE 
         CLASS(FTStack)                     :: self
         CLASS(FTObject)          , POINTER :: obj
         CLASS(FTLinkedListRecord), POINTER :: newRecord => NULL()
         CLASS(FTLinkedListRecord), POINTER :: tmp       => NULL()
         
         ALLOCATE(newRecord)
         CALL newRecord % initWithObject(obj)
         
         IF ( .NOT.ASSOCIATED(self % head) )     THEN
            self % head => newRecord
            self % tail => newRecord
         ELSE
            tmp                => self % head
            self % head        => newRecord
            self % head % next => tmp
            tmp  % previous    => newRecord
         END IF
         self % nRecords = self % nRecords + 1
         
      END SUBROUTINE push
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION peek(self)
!
!        -----------------------------------------
!        Return the object at the head of the list
!        ** No change of ownership **
!        -----------------------------------------
!
         IMPLICIT NONE 
         CLASS(FTStack)           :: self
         CLASS(FTObject), POINTER :: peek
         
         IF ( .NOT. ASSOCIATED(self % head) )     THEN
            peek => NULL()
            RETURN 
         END IF 

         peek => self % head % recordObject

      END FUNCTION peek    
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE pop(self,p)
!
!        ---------------------------------------------------
!        Remove the head of the list and return the object
!        that it points to. Calling routine gains ownership 
!        of the object.
!        ---------------------------------------------------
!
         IMPLICIT NONE  
         CLASS(FTStack)                     :: self
         CLASS(FTObject)          , POINTER :: p
         CLASS(FTLinkedListRecord), POINTER :: tmp => NULL()
         
         IF ( .NOT. ASSOCIATED(self % head) )     THEN
            p => NULL()
            RETURN 
         END IF 
            
         p => self % head % recordObject
         IF(.NOT.ASSOCIATED(p)) RETURN 
         CALL p % retain()
         
         tmp => self % head
         self % head => self % head % next
         
         CALL release(tmp)
         self % nRecords = self % nRecords - 1

      END SUBROUTINE pop
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION stackFromObject(obj) RESULT(cast)
!
!     -----------------------------------------------------
!     Cast the base class FTObject to the LinkedList class
!     -----------------------------------------------------
!
         IMPLICIT NONE  
         CLASS(FTObject), POINTER :: obj
         CLASS(FTStack) , POINTER :: cast
         
         cast => NULL()
         SELECT TYPE (e => obj)
            TYPE is (FTStack)
               cast => e
            CLASS DEFAULT
               
         END SELECT
         
      END FUNCTION stackFromObject
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE printStackDescription(self, iUnit) 
         IMPLICIT NONE 
         CLASS(FTStack) :: self
         INTEGER        :: iUnit
         
         CALL self % FTLinkedList % printDescription(iUnit = iUnit)
         
      END SUBROUTINE printStackDescription
!
!//////////////////////////////////////////////////////////////////////// 
! 
!      -----------------------------------------------------------------
!> Class name returns a string with the name of the type of the object
!>
!>  ### Usage:
!>
!>        PRINT *,  obj % className()
!>        if( obj % className = "FTStack")
!>
      FUNCTION stackClassName(self)  RESULT(s)
         IMPLICIT NONE  
         CLASS(FTStack)                             :: self
         CHARACTER(LEN=CLASS_NAME_CHARACTER_LENGTH) :: s
         
         s = "FTStack"
 
      END FUNCTION stackClassName
    
      END Module FTStackClass    