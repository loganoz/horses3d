!
!////////////////////////////////////////////////////////////////////////
!
!      FTMutableObjectArray.f90
!      Created: February 7, 2013 3:24 PM 
!      By: David Kopriva  
!
!>FTMutableObjectArray is a mutable array class to which objects
!>can be added, removed, replaced and accessed according to their 
!>index in the array.
!>
!>Fortran has pointers to arrays, but not arrays of pointers. To do the latter, one creates
!>a wrapper derived type and creates an array of that wrapper type. Fortran arrays are great, but
!>they are of fixed length, and they don't easily implement reference counting to keep track of
!>memory. For that, we have the FTMutableObjectArray. Performance reasons dictate that you 
!>will use regular arrays for numeric types and the like, but for generic objects we would use
!>an Object Array.
!>
!>You initialize a FTMutableObjectArray with the number of objects that you expect it to hold.
!>However, it can re-size itself if necessary. To be efficient, it adds more than one entry at a time
!>given by the ``chunkSize'', which you can choose for yourself. (The default is 10.)
!>##Definition
!>           TYPE(FTMutableObjectArray) :: array
!>#Usage
!>##Initialization
!>      CLASS(FTMutableObjectArray)  :: array
!>      INTEGER                      :: N = 11
!>      CALL array % initWithSize(N)
!>#Destruction
!>           CALL array  %  destuct() [Non Pointers]
!>           call release(array) [Pointers]
!>#Adding an Object
!>           TYPE(FTObject) :: obj
!>           obj => r1
!>           CALL array % addObject(obj)
!>#Removing an Object
!>           TYPE(FTObject) :: obj
!>           CALL array % removeObjectAtIndex(i)
!>#Accessing an Object
!>           TYPE(FTObject) :: obj
!>           obj => array % objectAtIndex(i)
!>#Replacing an Object
!>           TYPE(FTObject) :: obj
!>           obj => r1
!>           CALL array % replaceObjectAtIndexWithObject(i,obj)
!>#Setting the Chunk Size
!>           CALL array % setChunkSize(size)
!>#Finding The Number Of Items In The Array
!>           n =  array % count()
!>#Finding The Actual Allocated Size Of The Array
!>           n =  array % allocatedSize()

!
!////////////////////////////////////////////////////////////////////////
!
      MODULE FTMutableObjectArrayClass
      USE FTObjectClass
      IMPLICIT NONE
      
      TYPE FTObjectPointerWrapper
         CLASS(FTObject), POINTER ::  object => NULL()
      END TYPE FTObjectPointerWrapper
      
      PRIVATE :: FTObjectPointerWrapper
      PRIVATE :: increaseArraysize
      
      TYPE, EXTENDS(FTObject) ::  FTMutableObjectArray
         INTEGER                                            , PRIVATE :: count_
         TYPE(FTObjectPointerWrapper), DIMENSION(:), POINTER, PRIVATE :: array => NULL()
         INTEGER                                            , PRIVATE :: chunkSize_ = 10
!
!        --------
         CONTAINS
!        --------
!
         PROCEDURE, PUBLIC :: initWithSize => initObjectArrayWithSize
         PROCEDURE, PUBLIC :: destruct     => destructObjectArray
         PROCEDURE, PUBLIC :: addObject    => addObjectToArray
         PROCEDURE, PUBLIC :: replaceObjectAtIndexWithObject
         PROCEDURE, PUBLIC :: removeObjectAtIndex
         PROCEDURE, PUBLIC :: objectAtIndex
         
         PROCEDURE, PUBLIC :: printDescription => printArray
         PROCEDURE, PUBLIC :: className        => arrayClassName
!         
         PROCEDURE, PUBLIC :: setChunkSize
         PROCEDURE, PUBLIC :: chunkSize
         PROCEDURE, PUBLIC :: COUNT => numberOfItems
         PROCEDURE, PUBLIC :: allocatedSize
         
      END TYPE 
         
      INTERFACE cast
         MODULE PROCEDURE castToMutableObjectArray
      END INTERFACE cast
      
      INTERFACE release
         MODULE PROCEDURE releaseFTMutableObjectArray 
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
!> the array remains empty.
!>
!> *Usage
!>
!>       CLASS(FTMutableObjectArray)  :: array
!>       integer                      :: N = 11
!>       CALL array % initWithSize(N)
!>
      SUBROUTINE initObjectArrayWithSize( self, arraySize )    
         IMPLICIT NONE  
         CLASS( FTMutableObjectArray) :: self
         INTEGER                      :: arraySize
         INTEGER                      :: i
         
         CALL self % FTObject % init()
         
         ALLOCATE( self % array(arraySize) )
         
         DO i = 1,  arraySize
            self % array(i) %  object => NULL()
         END DO
         self % count_ = 0
         
      END SUBROUTINE initObjectArrayWithSize
!
!//////////////////////////////////////////////////////////////////////// 
! 
!>
!> Destructor for the class. This is called automatically when the
!> reference count reaches zero. Do not call this yourself.
!>
      SUBROUTINE destructObjectArray(self)  
         IMPLICIT NONE
         CLASS( FTMutableObjectArray) :: self
         CLASS(FTObject), POINTER     :: obj     => NULL()
         INTEGER                      :: i

         DO i = 1, self % count_
            obj => self % array(i) % object 
            IF ( ASSOCIATED(obj) ) CALL releaseFTObject(self = obj)
         END DO
         
         DEALLOCATE(self % array)
         self % array => NULL()
         self % count_ = 0  

      END SUBROUTINE
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
      SUBROUTINE releaseFTMutableObjectArray(self)  
         IMPLICIT NONE
         TYPE(FTMutableObjectArray) , POINTER :: self
         CLASS(FTObject), POINTER :: obj
         
         IF(.NOT. ASSOCIATED(self)) RETURN
         
         obj => self
         CALL releaseFTObject(self = obj)
         IF ( .NOT. ASSOCIATED(obj) )     THEN
            self => NULL() 
         END IF      
      END SUBROUTINE releaseFTMutableObjectArray
!
!//////////////////////////////////////////////////////////////////////// 
! 
!>
!> Add an object to the end of the array
!>
!> *Usage
!>
!>       CLASS(FTMutableObjectArray)      :: array
!>       CLASS(FTObject)        , POINTER :: obj
!>       CLASS(FTObjectSubclass), POINTER :: p
!>       obj => p
!>       CALL array % addObject(obj)
!>
      SUBROUTINE addObjectToArray(self,obj)
         IMPLICIT NONE  
         CLASS(FTMutableObjectArray) :: self
         CLASS(FTObject), POINTER    :: obj

         self % count_ = self % count_ + 1

         IF ( self % count_ > SIZE(self % array) )     THEN
            CALL increaseArraysize( self, self % count_ ) 
         END IF 
         
         self % array(self % count_) %  object => obj
         CALL obj % retain()
         
      END SUBROUTINE addObjectToArray
!
!//////////////////////////////////////////////////////////////////////// 
! 
!>
!> Remove an object at the index indx
!>
!> *Usage
!>
!>       CLASS(FTMutableObjectArray) :: array
!>       INTEGER                     :: indx
!>       CALL array % removeObjectAtIndex(indx)
!>
      SUBROUTINE removeObjectAtIndex(self,indx)  
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(FTMutableObjectArray) :: self
         INTEGER                     :: indx
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                     :: i
         CLASS(FTObject), POINTER    :: obj  => NULL()
         
         obj => self % array(indx) %  object
         
         IF ( ASSOCIATED(obj) )     THEN
            CALL releaseFTObject(self = obj)
         END IF 
         
         DO i = indx, self % count_-1
            self % array(i) % object => self % array(i+1) % object
         END DO
         self % array(self % count_) % object => NULL()
         self % count_                        = self % count_ - 1
         
      END SUBROUTINE removeObjectAtIndex
!
!//////////////////////////////////////////////////////////////////////// 
! 
!>
!> Replace an object at the index indx
!>
!> Usage
!> -----
!>
!>       CLASS(FTMutableObjectArray) :: array
!>       INTEGER                     :: indx
!>       CALL array % replaceObjectAtIndexWithObject(indx)
!>
      SUBROUTINE replaceObjectAtIndexWithObject(self,indx,replacement)  
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(FTMutableObjectArray) :: self
         INTEGER                     :: indx
         CLASS(FTObject), POINTER    :: replacement
!
!        ---------------
!        Local variables
!        ---------------
!
         CLASS(FTObject), POINTER    :: obj => NULL()
         
         obj => self % array(indx) %  object
         
         CALL releaseFTObject(obj)
         
         self % array(indx) %  object => replacement
         CALL replacement % retain()
         
      END SUBROUTINE  replaceObjectAtIndexWithObject
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE printArray(self,iUnit)  
         IMPLICIT NONE  
         CLASS(FTMutableObjectArray) :: self
         INTEGER                     :: iUnit
         INTEGER                     :: i
         CLASS(FTObject), POINTER    :: obj => NULL()
         
         DO i = 1, self % count_
            obj => self % array(i) % object
            CALL obj % printDescription(iUnit)
         END DO  
      END SUBROUTINE printArray
!
!//////////////////////////////////////////////////////////////////////// 
! 
!>
!> Access the object at the index indx
!>
!> *Usage
!>
!>       CLASS(FTMutableObjectArray) :: array
!>       INTEGER                     :: indx
!>       CLASS(FTObject), POINTER    :: obj
!>       obj => array % objectAtIndex(indx)
!>
      FUNCTION objectAtIndex(self,indx)  RESULT(obj)
         IMPLICIT NONE  
         CLASS(FTMutableObjectArray) :: self
         INTEGER                     :: indx
         CLASS(FTObject), POINTER    :: obj
         IF ( indx > self % count_ .OR. indx < 1)     THEN
            obj => NULL() 
         ELSE 
            obj => self % array(indx) % object
         END IF 
         
      END FUNCTION objectAtIndex
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE increaseArraySize( self, n )  
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         CLASS( FTMutableObjectArray) :: self
         INTEGER                      :: n
!
!        ---------------
!        Local Variables
!        ---------------
!
         TYPE(FTObjectPointerWrapper), DIMENSION(:), POINTER :: newArray 
         INTEGER                                             :: i, m
         
         IF ( n <= SIZE(self % array) )     THEN
            RETURN 
         END IF 
         
         m = (n - SIZE(self % array))/self % chunkSize_ + 1
         ALLOCATE( newArray(SIZE(self % array) + m*self % chunkSize_) )
         
         DO i = 1,  SIZE(self % array)
            newArray(i) %  object => self % array(i) %  object
         END DO
         DO i = SIZE(self % array) + 1,  SIZE(newArray)
            newArray(i) %  object => NULL()
         END DO
         
         DEALLOCATE(self % array)
         self % array => newArray

      END SUBROUTINE increaseArraySize
!
!//////////////////////////////////////////////////////////////////////// 
! 
!>
!> Set the number of items to be added when the array needs to be re-sized
!>
!> *Usage
!>
!>       CLASS(FTMutableObjectArray) :: array
!>       INTEGER                     :: sze = 42
!>       CALL array % setChunkSize(sze)
!>
      SUBROUTINE setChunkSize(self,chunkSize)  
         IMPLICIT NONE  
         CLASS( FTMutableObjectArray) :: self
         INTEGER              :: chunkSize
         
         self % chunkSize_ = chunkSize
         
      END SUBROUTINE setChunkSize
!
!//////////////////////////////////////////////////////////////////////// 
! 
!>
!> Returns the number of items to be added when the array needs to be re-sized
!>
!> *Usage
!>
!>       CLASS(FTMutableObjectArray) :: array
!>       INTEGER                     :: sze
!>       sze =  array % chunkSize
!>
      INTEGER FUNCTION chunkSize(self)  
         IMPLICIT NONE  
         CLASS( FTMutableObjectArray) :: self
         chunkSize = self % chunkSize_
      END FUNCTION chunkSize
!
!//////////////////////////////////////////////////////////////////////// 
! 
!>
!> Generic name: count
!>
!> Returns the actual number of items in the array.
!>
!> *Usage
!>
!>       CLASS(FTMutableObjectArray) :: array
!>       INTEGER                     :: sze
!>       sze =  array % count()
!>
      INTEGER FUNCTION numberOfItems(self)  
         IMPLICIT NONE  
         CLASS( FTMutableObjectArray) :: self
         numberOfItems = self % count_
      END FUNCTION numberOfItems
!
!//////////////////////////////////////////////////////////////////////// 
! 
      INTEGER FUNCTION allocatedSize(self)  
         IMPLICIT NONE  
         CLASS( FTMutableObjectArray) :: self
         IF ( ASSOCIATED(self % array) )     THEN
            allocatedSize = SIZE(self % array)
         ELSE
            allocatedSize = 0
         END IF 
      END FUNCTION allocatedSize
!
!---------------------------------------------------------------------------
!> Generic Name: cast
!> 
!> Cast a pointer to the base class to an FTMutableObjectArray pointer 
!---------------------------------------------------------------------------
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION objectArrayFromObject(obj) RESULT(cast)
!
!     -----------------------------------------------------
!     Cast the base class FTObject to the FTException class
!     -----------------------------------------------------
!
         IMPLICIT NONE  
         CLASS(FTObject)            , POINTER :: obj
         CLASS(FTMutableObjectArray), POINTER :: cast
         
         cast => NULL()
         SELECT TYPE (e => obj)
            TYPE is (FTMutableObjectArray)
               cast => e
            CLASS DEFAULT
               
         END SELECT
         
      END FUNCTION objectArrayFromObject
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE castToMutableObjectArray(obj,cast) 
!
!     -----------------------------------------------------
!     Cast the base class FTObject to the FTException class
!     -----------------------------------------------------
!
         IMPLICIT NONE  
         CLASS(FTObject)            , POINTER :: obj
         CLASS(FTMutableObjectArray), POINTER :: cast
         
         cast => NULL()
         SELECT TYPE (e => obj)
            TYPE is (FTMutableObjectArray)
               cast => e
            CLASS DEFAULT
               
         END SELECT
         
      END SUBROUTINE castToMutableObjectArray
!
!//////////////////////////////////////////////////////////////////////// 
! 
!      -----------------------------------------------------------------
!> Class name returns a string with the name of the type of the object
!>
!>  ### Usage:
!>
!>        PRINT *,  obj % className()
!>        if( obj % className = "FTMutableObjectArray")
!>
      FUNCTION arrayClassName(self)  RESULT(s)
         IMPLICIT NONE  
         CLASS(FTMutableObjectArray)                :: self
         CHARACTER(LEN=CLASS_NAME_CHARACTER_LENGTH) :: s
         
         s = "FTMutableObjectArray"
 
      END FUNCTION arrayClassName

      
      END Module  FTMutableObjectArrayClass    