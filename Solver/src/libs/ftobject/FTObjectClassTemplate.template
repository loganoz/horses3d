!
!////////////////////////////////////////////////////////////////////////
!
!      <ClassName>.f90
!
!><ClassName> is a ...
!>
!>##Definition
!>           TYPE(<ClassName>) :: varName
!>#Usage
!>##Initialization
!>      CLASS(<ClassName>)  :: FTMutable<ClassName>
!>      INTEGER                      :: N= 11
!>      CALL FTMutable<ClassName> % init...
!>#Destruction
!>           CALL FTMutable<ClassName>  %  destuct() [Non Pointers]
!>           call release(FTMutable<ClassName>) [Pointers]

!
!////////////////////////////////////////////////////////////////////////
!
      MODULE <ClassName>Class
      USE FTObjectClass
      IMPLICIT NONE
      
      TYPE, EXTENDS(FTObject) ::  <ClassName>
      !Instance variables
!
!        --------
         CONTAINS
!        --------
!
         PROCEDURE, PUBLIC :: init<ClassName> => init<ClassName><ClassName>
         PROCEDURE, PUBLIC :: destruct     => destruct<ClassName>
         
         PROCEDURE, PUBLIC :: printDescription => print<ClassName>
         PROCEDURE, PUBLIC :: className        => <ClassName>ClassName
!         
         PROCEDURE, PUBLIC :: setChunkSize
         PROCEDURE, PUBLIC :: chunkSize
         PROCEDURE, PUBLIC :: COUNT => numberOfItems
         PROCEDURE, PUBLIC :: allocatedSize
         
      END TYPE 
         
      INTERFACE cast
         MODULE PROCEDURE castTo<ClassName>
      END INTERFACE cast
      
      INTERFACE release
         MODULE PROCEDURE release<ClassName> 
      END INTERFACE  
!
!     ======== 
      CONTAINS  
!     ======== 
!
!
!//////////////////////////////////////////////////////////////////////// 
! 
!>
!> Designated initializer. Initializes the amount of storage, but
!> the <ClassName> remains empty.
!>
!> *Usage
!>
!>       CLASS(<ClassName>)  :: <ClassName>
!>       integer                      :: N = 11
!>       CALL <ClassName> % init<ClassName>(N)
!>
      SUBROUTINE init<ClassName><ClassName>( self, <ClassName>Size )    
         IMPLICIT NONE  
         CLASS( <ClassName>) :: self
         INTEGER                      :: <ClassName>Size
         INTEGER                      :: i
         
         CALL self % FTObject % init()
         
         
      END SUBROUTINE init<ClassName><ClassName>
!
!//////////////////////////////////////////////////////////////////////// 
! 
!>
!> Destructor for the class. This is called automatically when the
!> reference count reaches zero. Do not call this yourself.
!>
      SUBROUTINE destruct<ClassName>(self)  
         IMPLICIT NONE
         CLASS( <ClassName>) :: self
         CLASS(FTObject), POINTER     :: obj     => NULL()

      END SUBROUTINE destruct<ClassName>
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
      SUBROUTINE release<ClassName>(self)  
         IMPLICIT NONE
         TYPE(<ClassName>) , POINTER :: self
         CLASS(FTObject), POINTER :: obj
         obj => self
         CALL releaseFTObject(self = obj)
         IF ( .NOT. ASSOCIATED(obj) )     THEN
            self => NULL() 
         END IF      
      END SUBROUTINE release<ClassName>
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE print<ClassName>(self,iUnit)  
         IMPLICIT NONE  
         CLASS(<ClassName>) :: self
         INTEGER                     :: iUnit
         INTEGER                     :: i
         CLASS(FTObject), POINTER    :: obj => NULL()
         
         DO i = 1, self % count_
            obj => self % <ClassName>(i) % object
            CALL obj % printDescription(iUnit)
         END DO  
      END SUBROUTINE print<ClassName>

!
!---------------------------------------------------------------------------
!> Generic Name: cast
!> 
!> Cast a pointer to the base class to an <ClassName> pointer 
!---------------------------------------------------------------------------
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION <ClassName>FromObject(obj) RESULT(cast)
!
!     -----------------------------------------------------
!     Cast the base class FTObject to the FTException class
!     -----------------------------------------------------
!
         IMPLICIT NONE  
         CLASS(FTObject)            , POINTER :: obj
         CLASS(<ClassName>), POINTER :: cast
         
         cast => NULL()
         SELECT TYPE (e => obj)
            TYPE is (<ClassName>)
               cast => e
            CLASS DEFAULT
               
         END SELECT
         
      END FUNCTION <ClassName>FromObject
!
!//////////////////////////////////////////////////////////////////////// 
! 
!      -----------------------------------------------------------------
!> Class name returns a string with the name of the type of the object
!>
!>  ### Usage:
!>
!>        PRINT *,  obj % className()
!>        if( obj % className = "<ClassName>")
!>
      FUNCTION <ClassName>ClassName(self)  RESULT(s)
         IMPLICIT NONE  
         CLASS(<ClassName>)                :: self
         CHARACTER(LEN=CLASS_NAME_CHARACTER_LENGTH) :: s
         
         s = "<ClassName>"
 
      END FUNCTION <ClassName>ClassName

      
      END Module  <ClassName>Class    