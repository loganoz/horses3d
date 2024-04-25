!
!////////////////////////////////////////////////////////////////////////
!
!      SparseMatrixClass.f90
!      Created: July 29, 2013 10:59 AM 
!      By: David Kopriva  
!
!
!////////////////////////////////////////////////////////////////////////
!
!>FTSparseMatrixData is used by the FTSparseMatrix Class. Users will 
!>usually not interact with or use this class directly.
!>
      Module FTSparseMatrixData 
      USE FTObjectClass
      IMPLICIT NONE
!
!     ---------------
!     Type definition
!     ---------------
!
      TYPE, EXTENDS(FTObject) :: MatrixData
         INTEGER                  :: key
         CLASS(FTObject), POINTER :: object
!
!        ========
         CONTAINS
!        ========
!
         PROCEDURE :: initWithObjectAndKey
         PROCEDURE :: destruct => destructMatrixData
         
      END TYPE MatrixData
      
      INTERFACE cast
         MODULE PROCEDURE castObjectToMatrixData
      END INTERFACE cast
      
!
!     ========      
      CONTAINS
!     ========
!
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE initWithObjectAndKey(self,object,key)
!
!        ----------------------
!        Designated initializer
!        ----------------------
!
         IMPLICIT NONE
         CLASS(MatrixData)        :: self
         CLASS(FTObject), POINTER :: object
         INTEGER                  :: key
         
         CALL self % FTObject % init()
         
         self % key = key
         self % object => object
         CALL self % object % retain()
         
      END SUBROUTINE initWithObjectAndKey
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE destructMatrixData(self)
         IMPLICIT NONE  
         CLASS(MatrixData) :: self
         
         IF ( ASSOCIATED(self % object) )     THEN
            CALL releaseFTObject(self = self % object)
         END IF 
         
         CALL self % FTObject % destruct

      END SUBROUTINE destructMatrixData
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE castObjectToMatrixData(obj,cast)  
         IMPLICIT NONE  
!
!        -----------------------------------------------------
!        Cast the base class FTObject to the FTException class
!        -----------------------------------------------------
!
         CLASS(FTObject)  , POINTER :: obj
         CLASS(MatrixData), POINTER :: cast
         
         cast => NULL()
         SELECT TYPE (e => obj)
            TYPE is (MatrixData)
               cast => e
            CLASS DEFAULT
               
         END SELECT
         
      END SUBROUTINE castObjectToMatrixData
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION matrixDataCast(obj)  RESULT(cast)
         IMPLICIT NONE  
!
!        -----------------------------------------------------
!        Cast the base class FTObject to the FTException class
!        -----------------------------------------------------
!
         CLASS(FTObject)  , POINTER :: obj
         CLASS(MatrixData), POINTER :: cast
         
         cast => NULL()
         SELECT TYPE (e => obj)
            TYPE is (MatrixData)
               cast => e
            CLASS DEFAULT
               
         END SELECT
         
      END FUNCTION matrixDataCast
      
      END Module FTSparseMatrixData
!@mark -
!>The sparse matrix stores an FTObject pointer associated
!>with two keys (i,j) as a hash table.
!>
!>Hash tables are data structures designed to enable storage and fast
!>retrieval of key-value pairs. An example of a key-value pair is
!>a variable name (``gamma'') and its associated value (``1.4'').
!>The table itself is typically an array.
!>The location of the value in a hash table associated with
!>a key, $k$, is specified by way of a hash function, $H(k)$.
!>In the case of a variable name and value, the hash function
!>would convert the name into an integer that tells where to
!>find the associated value in the table.
!>
!>A very simple example of a
!>hash table is, in fact, a singly dimensioned array. The key is 
!>the array index and the value is what is stored at that index.
!>Multiple keys can be used to identify data; a two dimensional
!>array provides an example of where two keys are used to access memory
!>and retrieve the value at that location.
!>If we view a singly dimensioned array as a special case of a hash table,
!>its hash function is just the array index, $H(j)=j$. A doubly dimensioned array
!>could be (and often is) stored columnwise as a singly dimensioned array by creating a hash
!>function that maps the two indices to a single location in the array, e.g.,
!>$H(i,j) = i + j*N$, where $N$ is the range of the first index, $i$. 
!>
!>Two classes are included in FTObjectLibrary. The first, FTSparseMatrix, works with an ordered pair, (i,j), as the
!>keys. The second, FTMultiIndexTable, uses an array of integers as the keys.
!>
!>Both classes include enquiry functions to see of an object exists for the given keys. Otherwise,
!>the function that returns an object for a given key will return an UNASSOCIATED pointer if there
!>is no object for the key. Be sure to retain any object returned by the objectForKeys methods if 
!>you want to keep it beyond the lifespan of the matrix or table. For example,
!>
!>           TYPE(FTObject) :: obj
!>           obj => matrix % objectForKeys(i,j)
!>           IF ( ASSOCIATED(OBJ) ) THEN
!>               CALL obj % retain()
!>                 Cast obj to something useful
!>           ELSE
!>              Perform some kind of error recovery
!>           END IF 
!>The sparse matrix stores an FTObject pointer associated
!>with two keys (i,j) as a hash table. The size, N = the range of i.
!>
!>##Definition (Subclass of FTObject)
!>
!>         TYPE(FTSparseMatrix) :: SparseMatrix
!>#Usage
!>##Initialization
!>
!>         CALL SparseMatrix % initWithSize(N)
!>
!>##Destruction
!>
!>         CALL release(SparseMatrix)
!>
!>##Adding an object
!>
!>         CLASS(FTObject), POINTER :: obj
!>         CALL SparseMatrix % addObjectForKeys(obj,i,j)
!>
!>##Retrieving an object
!>
!>         CLASS(FTObject), POINTER :: obj
!>         obj => SparseMatrix % objectForKeys(i,j)
!>
!>Be sure to retain the object if you want it to live
!>      beyond the life of the table.
!>
!>##Testing the presence of keys
!>
!>         LOGICAL :: exists
!>         exists = SparseMatrix % containsKeys(i,j)
!
!////////////////////////////////////////////////////////////////////////
!
      Module FTSparseMatrixClass
      USE FTObjectClass
      USE FTLinkedListClass
      USE FTLinkedListIteratorClass
      USE FTSparseMatrixData
      IMPLICIT NONE
!
!     ----------------------
!     Class type definitions
!     ----------------------
!
      TYPE FTLinkedListPtr
         CLASS(FTLinkedList), POINTER :: list
      END TYPE FTLinkedListPtr
      PRIVATE :: FTLinkedListPtr
      
      TYPE, EXTENDS(FTObject) :: FTSparseMatrix
         TYPE(FTLinkedListPtr)     , DIMENSION(:), ALLOCATABLE :: table
         TYPE(FTLinkedListIterator), PRIVATE                   :: iterator
!
!        ========
         CONTAINS
!        ========
!
         PROCEDURE :: initWithSize     => initSparseMatrixWithSize
         PROCEDURE :: destruct         => destructSparseMatrix
         PROCEDURE :: containsKeys     => SparseMatrixContainsKeys
         PROCEDURE :: addObjectForKeys => addObjectToSparseMatrixForKeys
         PROCEDURE :: objectForKeys    => objectInSparseMatrixForKeys
         PROCEDURE :: SparseMatrixSize
         
      END TYPE FTSparseMatrix
      
      INTERFACE release
         MODULE PROCEDURE releaseFTSparseMatrix 
      END INTERFACE  
      
!
!     ========
      CONTAINS
!     ========
!
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE initSparseMatrixWithSize(self,N)  
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(FTSparseMatrix) :: self
         INTEGER            :: N
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: j
         
         CALL self % FTObject % init()
         
         ALLOCATE(self % table(N))
         DO j = 1, N
            ALLOCATE(self % table(j) % list)
            CALL self % table(j) % list % init()
         END DO
         
         CALL self % iterator % init()
         
      END SUBROUTINE initSparseMatrixWithSize
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE addObjectToSparseMatrixForKeys(self,obj,i,j)
         IMPLICIT NONE  
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(FTSparseMatrix)    :: self
         CLASS(FTObject), POINTER :: obj
!
!        ---------------
!        Local variables
!        ---------------
!
         CLASS(MatrixData), POINTER :: mData
         CLASS(FTObject)  , POINTER :: ptr
         INTEGER                    :: i,j
         
         IF ( .NOT.self % containsKeys(i,j) )     THEN
            ALLOCATE(mData)
            CALL mData % initWithObjectAndKey(obj,j)
            ptr => mData
            CALL self % table(i) % list % add(ptr)
            CALL releaseFTObject(ptr)
         END IF 
         
      END SUBROUTINE addObjectToSparseMatrixForKeys
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION objectInSparseMatrixForKeys(self,i,j) RESULT(r)
!
!     ---------------------------------------------------------------
!     Returns the stored FTObject for the keys (i,j). Returns NULL()
!     if the object isn't in the table. Retain the object if it needs
!     a strong reference by the caller.
!     ---------------------------------------------------------------
!
         IMPLICIT NONE  
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(FTSparseMatrix)    :: self
         INTEGER                  :: i,j
         CLASS(FTObject), POINTER :: r
!
!        ---------------
!        Local variables
!        ---------------
!
         CLASS(MatrixData)  , POINTER :: mData
         CLASS(FTObject)    , POINTER :: obj
         CLASS(FTLinkedList), POINTER :: list
         
         r    => NULL()
         IF(.NOT.ALLOCATED(self % table))     RETURN 
         list => self % table(i) % list
         IF(.NOT.ASSOCIATED(list))    RETURN 
         IF (  list % COUNT() == 0 )  RETURN
!
!        ----------------------------
!        Step through the linked list
!        ----------------------------
!
         r => NULL()
         
         CALL self % iterator % setLinkedList(self % table(i) % list)
         DO WHILE (.NOT.self % iterator % isAtEnd())
         
            obj => self % iterator % object()
            CALL cast(obj,mData)
            IF ( mData % key == j )     THEN
               r => mData % object
               EXIT 
            END IF 
            
            CALL self % iterator % moveToNext()
         END DO

      END FUNCTION objectInSparseMatrixForKeys
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION SparseMatrixContainsKeys(self,i,j)  RESULT(r)
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(FTSparseMatrix) :: self
         INTEGER                :: i, j
         LOGICAL                :: r
!
!        ---------------
!        Local variables
!        ---------------
!
         CLASS(FTObject)    , POINTER :: obj
         CLASS(MatrixData)  , POINTER :: mData
         CLASS(FTLinkedList), POINTER :: list
         
         r = .FALSE.
         IF(.NOT.ALLOCATED(self % table))                RETURN 
         IF(.NOT.ASSOCIATED(self % table(i) % list))     RETURN
         IF ( self % table(i) % list % COUNT() == 0 )    RETURN 
!
!        ----------------------------
!        Step through the linked list
!        ----------------------------
!
         list => self % table(i) % list
         CALL self % iterator % setLinkedList(list)
         CALL self % iterator % setToStart()
         DO WHILE (.NOT.self % iterator % isAtEnd())
         
            obj => self % iterator % object()
            CALL cast(obj,mData)
            IF ( mData % key == j )     THEN
               r = .TRUE.
               RETURN  
            END IF 
            
            CALL self % iterator % moveToNext()
         END DO
         
      END FUNCTION SparseMatrixContainsKeys
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE destructSparseMatrix(self)
         IMPLICIT NONE  
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(FTSparseMatrix) :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: j
         
         DO j = 1, SIZE(self % table)
            IF ( ASSOCIATED(self % table(j) % list) )     THEN
               CALL releaseSMMemberList(list = self % table(j) % list)
            END IF 
         END DO

         IF(ALLOCATED(self % table))   DEALLOCATE(self % table)

         CALL self % iterator % destruct()
         
         CALL self % FTObject % destruct()
         
      END SUBROUTINE destructSparseMatrix
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
      SUBROUTINE releaseFTSparseMatrix(self)  
         IMPLICIT NONE
         CLASS(FTSparseMatrix) , POINTER :: self
         CLASS(FTObject)       , POINTER :: obj
         
         IF(.NOT. ASSOCIATED(self)) RETURN
         
         obj => self
         CALL releaseFTObject(self = obj)
         IF ( .NOT. ASSOCIATED(obj) )     THEN
            self => NULL() 
         END IF      
      END SUBROUTINE releaseFTSparseMatrix
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE releaseSMMemberList(list)  
          IMPLICIT NONE  
          CLASS(FTLinkedList), POINTER :: list
          CLASS(FTObject)    , POINTER :: obj

          obj => list
          CALL releaseFTObject(self = obj)
          IF(.NOT. ASSOCIATED(obj)) list => NULL()
      END SUBROUTINE releaseSMMemberList
!
!//////////////////////////////////////////////////////////////////////// 
! 
      INTEGER FUNCTION SparseMatrixSize(self)  
         IMPLICIT NONE  
         CLASS(FTSparseMatrix) :: self
         IF ( ALLOCATED(self % table) )     THEN
            SparseMatrixSize =  SIZE(self % table)
         ELSE
            SparseMatrixSize = 0
         END IF 
      END FUNCTION SparseMatrixSize
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION SparseMatrixFromObject(obj) RESULT(cast)
!
!     -----------------------------------------------------
!     Cast the base class FTObject to the FTException class
!     -----------------------------------------------------
!
         IMPLICIT NONE  
         CLASS(FTObject)   , POINTER :: obj
         CLASS(FTSparseMatrix), POINTER :: cast
         
         cast => NULL()
         SELECT TYPE (e => obj)
            TYPE is (FTSparseMatrix)
               cast => e
            CLASS DEFAULT
               
         END SELECT
         
      END FUNCTION SparseMatrixFromObject
!
!////////////////////////////////////////////////////////////////////////
!
      INTEGER FUNCTION Hash1( idPair )
         INTEGER, DIMENSION(2) :: idPair
         Hash1 = MAXVAL(idPair)
      END FUNCTION Hash1
!
!////////////////////////////////////////////////////////////////////////
!
      INTEGER FUNCTION Hash2( idPair )
         INTEGER, DIMENSION(2) :: idPair
         Hash2 = MINVAL(idPair)
      END FUNCTION Hash2

      END Module FTSparseMatrixClass