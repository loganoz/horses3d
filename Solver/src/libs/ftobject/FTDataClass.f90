!
!////////////////////////////////////////////////////////////////////////
!
!      FTDataClass.f90
!      Created: July 11, 2013 2:00 PM 
!      By: David Kopriva  
!
!>FTData defines a subclass of FTObject to contain immutable
!>generic data, including derived types. 
!>
!>The initializer
!>copies the data and takes ownership of that copy. FTData
!>gives a way to use derived types without having to subclass
!>FTObject.
!
!////////////////////////////////////////////////////////////////////////
!
      Module FTDataClass 
      USE FTObjectClass
      IMPLICIT NONE
!
!     ---------
!     Constants
!     ---------
!
      INTEGER, PARAMETER          :: DATA_CLASS_TYPE_LENGTH = 32
!
!     ---------------------
!     Class type definition
!     ---------------------
!
      TYPE, EXTENDS(FTObject) :: FTData
         PRIVATE 
         CHARACTER(LEN=DATA_CLASS_TYPE_LENGTH) :: dataType
         CHARACTER(LEN=1), ALLOCATABLE         :: dataStorage(:) 
!
!        ========         
         CONTAINS 
!        ========
!         
         PROCEDURE, PUBLIC :: initWithDataOfType
         PROCEDURE, PUBLIC :: storedData
         PROCEDURE, PUBLIC :: className => dataClassName
      END TYPE FTData
      
      INTERFACE release
         MODULE PROCEDURE releaseFTData 
      END INTERFACE  
      
      CONTAINS 
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE initWithDataOfType(self,genericData,dataType)  
         IMPLICIT NONE  
         CLASS(FTData)    :: self
         CHARACTER(LEN=*) :: dataType
         CHARACTER(LEN=1) :: genericData(:)
         
         INTEGER          :: dataSize
          
          CALL self % FTObject % init()
          
          dataSize = SIZE(genericData)
          ALLOCATE(self % dataStorage(dataSize))
          
          self % dataStorage = genericData
          self % dataType    = dataType
          
      END SUBROUTINE initWithDataOfType
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE releaseFTData(self)  
         IMPLICIT NONE
         CLASS(FTData)  , POINTER :: self
         CLASS(FTObject), POINTER :: obj
         
         IF(.NOT. ASSOCIATED(self)) RETURN
         
         obj => self
         CALL releaseFTObject(self = obj)
         IF ( .NOT. ASSOCIATED(obj) )     THEN
            self => NULL() 
         END IF      
      END SUBROUTINE releaseFTData
!@mark -
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION storedData(self)  RESULT(d)
         IMPLICIT NONE  
         CLASS(FTData)    :: self
         CHARACTER(LEN=1) :: d(SIZE(self%dataStorage))
         d = self % dataStorage
      END FUNCTION storedData
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION dataType(self)  RESULT(t)
         IMPLICIT NONE  
         CLASS(FTData)    :: self
         CHARACTER(LEN=DATA_CLASS_TYPE_LENGTH) :: t
         t = self % dataType
      END FUNCTION dataType
!
!//////////////////////////////////////////////////////////////////////// 
! 
!      -----------------------------------------------------------------
!> Class name returns a string with the name of the type of the object
!>
!>  ### Usage:
!>
!>        PRINT *,  obj % className()
!>        if( obj % className = "FTData")
!>
      FUNCTION dataClassName(self)  RESULT(s)
         IMPLICIT NONE  
         CLASS(FTData)                              :: self
         CHARACTER(LEN=CLASS_NAME_CHARACTER_LENGTH) :: s
         
         s = "FTData"
 
      END FUNCTION dataClassName
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION dataIsOfType(self, dataType)  RESULT(t)
         IMPLICIT NONE  
         CLASS(FTData)                         :: self
         CHARACTER(LEN=DATA_CLASS_TYPE_LENGTH) :: dataType
         LOGICAL                               :: t
         
         IF ( dataType == self % dataType )     THEN
            t = .TRUE. 
         ELSE 
            t = .FALSE. 
         END IF 
      END FUNCTION dataIsOfType
      
      END Module FTDataClass