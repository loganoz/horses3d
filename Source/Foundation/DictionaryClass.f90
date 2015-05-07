!
! /////////////////////////////////////////////////////////////////////
!
!
!     DictionaryClass
!
!!     Created on: June 14, 2002
!!
!!     Modification History:
!!        June 18, 2002 : Add reSize and error handler
!!        Nov 26, 2002  : Use a simple Hash table for lookup
!!        June 30, 2011 : Simplified version without case sensitivity or error handler.
!!
!!     This class defines a mutable dictionary. To the user it defines
!      a keyword/value relationship. 
!!   
!
!      TYPE Dictionary
!
!
!      PUBLIC METHODS:
!
!         Constructor/Destructor:
!
!            SUBROUTINE ConstructDictionary( this )
!            SUBROUTINE ConstructDictionary( this, type )
!            SUBROUTINE DestructDictionary ( this )
!
!         Adding entries to dictionary:
!
!            SUBROUTINE AddValue_ForKey_ToDict_( value, key, this )
!
!         Getting values from dictionary:
!
!            SUBROUTINE GetValue_ForKey_FromDict_( val, key, dictionary )
!               val can be INTEGER, REAL( RP ), CHARACTER or LOGICAL
!
!         Debugging
!
!            SUBROUTINE PrintDictionary( this, unit )
!         
!         Inquiry:
!
!            LOGICAL FUNCTION DictionaryIsConstructed( this )
!            LOGICAL FUNCTION KeyFoundInDictionary( key, dictionary )
!            INTEGER FUNCTION NumberOfEntriesInDictionary( dictionary )
!            TYPE(KeyValuePair) FUNCTION KeyValuePairForKeyFromDict( key, this )
!
!         Iterator Methods:
!
!            SUBROUTINE InitDictionaryIterator(this)
!            TYPE(KeyValuePair), POINTER FUNCTION DictionaryIteratorNext(this)
!
!      PRIVATE METHODS:
!
!         INTEGER FUNCTION Hash(s,n)
!                     
!!
!!     Usage:
!!        To create a dictionary, call the constructor. When it is
!!        no longer needed, call the destructor. To add keyword pairs, 
!!        call the AddValue method. Note that, currently, the value
!!        is assumed to stored as a string of length DICT_VALUE_STRING_LENGTH, 
!!        so if computed values are to be stored in a dictionary, they
!!        need to be converted to a string first before being added. Finally, 
!!        to access a dictionary entry, use GetValue_ForKey_FromDict_. One can
!!        Iterate through a dictionary to find all keys/values using the
!!        iterator methods.
!!
!!     EXAMPLE:
!!
!!>    TYPE( Dictionary )                        :: theDict
!!     CHARACTER( LEN=DICT_KWD_STRING_LENGTH )   :: kw
!!     DOUBLE PRECISION                          :: x
!!     TYPE(KeyValuePair), POINTER               :: p
!!
!!     CALL ConstructDictionary( theDict )
!!     
!!     kw = "aKeyword"
!!     x  = 3.14
!!     CALL AddValue_ForKey_ToDict_( x, kw, theDict )
!!
!!     CALL GetValue_ForKey_FromDict_( x, kw, theDict )
!!     WRITE( 6, * ) x
!!
!!     CALL InitDictionaryIterator(theDict)
!!     DO
!!        p => DictionaryIteratorNext(theDict)
!!        IF ( .NOT.ASSOCIATED(p) )     EXIT
!!        -- DO Something with p --
!!     END DO
!!<
!
!
!!    @author David A. Kopriva
!
! /////////////////////////////////////////////////////////////////////
!
!
!     DictionaryClass.F
!
!!
!!     Modification History:
!!        version 0.0 Nov 27, 2002
!!
!     This file contains two class definitions and one type difinition.
!     the fundamental type within which the entries are stored is the
!     KeyValuePair. The first class is an array of KeyValuePair's. It
!     should not be needed by uses of the dictionary class. The second
!     class is the dictionary itself. See below for documentation on how
!     to use a dictionary.
!
!     Classes Defined Here:
!
!        KeyValueArray
!        DictionaryClass
!
!     Types defined here
!
!      TYPE KeyValuePair
!      TYPE KeyValueArray
!      TYPE Dictionary
!
!      PUBLIC DATA:
!
!         INTEGER, PARAMETER, PUBLIC :: DICT_VALUE_STRING_LENGTH
!         INTEGER, PARAMETER, PUBLIC :: DICT_KWD_STRING_LENGTH
!         INTEGER, PARAMETER, PUBLIC :: CASE_SENSITIVE_DICTIONARY
!         INTEGER, PARAMETER, PUBLIC :: CASE_INSENSITIVE_DICTIONARY
!
!!    @author David A. Kopriva
!!    
! /////////////////////////////////////////////////////////////////////
!
!! The KeyValueArrayClass sefines a KeyValueArray which is used as
!! entries in a dictionary. It should not be needed outside of the
!! dictionary class.
!
!  Methods:
!
!     Construction/Destruction:
!
!        SUBROUTINE ConstructKeyValueArray( this, size )
!        SUBROUTINE DestructKeyValueArray( this )
!        LOGICAL FUNCTION KeyValueArrayIsConstructed(this)
!
!     Adding an item:
!        
!        SUBROUTINE AddItemToKeyValueArray( item, this )
!
!     Re-Sizing the array:
!        
!        SUBROUTINE ReSize( this )   (PRIVATE)
!
!     Accessing items:
!
!        TYPE(KeyValuePair) FUNCTION ItemForKeyInList( key, this ) RESULT(pair)

!     Iterator Methods
!
!        SUBROUTINE InitKeyValueArrayIterator(this)
!        FUNCTION KeyValueArrayIteratorNext(this) RESULT(p)
!        SUBROUTINE PrintKeyValueArray( this, unit )
!    
!!
!  ******
   MODULE KeyValueArrayClass
!  ******
!
     USE SMConstants
     IMPLICIT NONE
     
     PRIVATE :: ReSize
!
!    -----------------
!    Module constants:
!    -----------------
!       -----------------------------------------------------------
!!      Length of strings used in a dictionary. VALUE refers to the
!!      number of characters to represent the value, and KWD refers
!!      to the number of characters allowed for a keyword.
!       -----------------------------------------------------------
!!
        INTEGER, PARAMETER, PUBLIC  :: DICT_VALUE_STRING_LENGTH = 132
        INTEGER, PARAMETER, PUBLIC  :: DICT_KWD_STRING_LENGTH   = 32
        INTEGER, PARAMETER, PRIVATE :: KV_ALLOC_SIZE            = 5
!
!    ----------------
!    Type definitions
!    ----------------
!
!       --------------------------------------------------------
!!      The KeyValuePair stores the key and the value as strings
!       --------------------------------------------------------
!
         TYPE KeyValuePair
            CHARACTER( LEN=DICT_KWD_STRING_LENGTH )   :: key
            CHARACTER( LEN=DICT_VALUE_STRING_LENGTH ) :: value
         END TYPE KeyValuePair
!
!       --------------------------------------------------------
!!      The KeyValueArray stores an array of KeyValuePairs. A
!!      Dictionary will have an array of these arrays.
!       --------------------------------------------------------
!
        TYPE KeyValueArray
           LOGICAL                                   :: constructed
           INTEGER                                   :: size
           INTEGER                                   :: noOfItems
           TYPE(KeyValuePair), DIMENSION(:), POINTER :: items
           INTEGER                                   :: iteratorItemID
        END TYPE KeyValueArray
!
!    ========
     CONTAINS
!    ========
!
!                           **************************
!                           Methods for KeyValuePair's
!                           **************************
!
!     /////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------------
!!    SetPair is simply an accessor method for a KeyValuePair. Given a key and
!!    a value, it sets the pair.
!     -------------------------------------------------------------------------
!
      SUBROUTINE SetPair( this, key, value )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(KeyValuePair) :: this
      CHARACTER(LEN=*)   :: key
      CHARACTER(LEN=*)   :: value
!
      this%key   = key
      this%value = value
      
      END SUBROUTINE SetPair
!
!     /////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------------
!!    Undefined checks to see if the key and value are set to "undefined" and
!!    returns true if they are.
!     -------------------------------------------------------------------------
!
      LOGICAL FUNCTION Undefined( this )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(KeyValuePair) :: this
!
      IF ( TRIM(this%key)    == "undefined"      .AND.                        &
     &     TRIM(this%value)  == "undefined" )     THEN
         Undefined = .TRUE.
      ELSE
         Undefined = .FALSE.
      END IF
      
      END FUNCTION Undefined
!
!
!                           ***************************
!                           Methods for KeyValuearray's
!                           ***************************
!
!     /////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------------
!!    Constructor for the KeyValueArray. The default constructor, without the
!!    optional size parameter will create a list of length KV_ALLOC_SIZE.
!     -------------------------------------------------------------------------
!
      SUBROUTINE ConstructKeyValueArray(this, size)
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(KeyValueArray) :: this
      INTEGER, OPTIONAL  :: size
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER                             :: s
      INTEGER                             :: listSize
      INTEGER                             :: j
!      
      this%constructed = .FALSE.
      IF ( PRESENT(size) )     THEN
         listSize = size
      ELSE
         listSize = KV_ALLOC_SIZE
      END IF
      
      ALLOCATE( this%items(listSize), STAT = s )
      
      IF ( s /= 0 )     THEN
         PRINT *, "Error Constructing KeyValueArray"
         RETURN
      END IF
      this%size        = listSize
      this%noOfItems   = 0
      this%constructed = .TRUE.
!
!     ----------------------------------
!     Initialize elements to "undefined"
!     ----------------------------------
!      
      DO j = 1, listSize
         CALL SetPair( this%items(j), "undefined", "undefined" )
      END DO
!
      END SUBROUTINE ConstructKeyValueArray
!
!     /////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------
!!    The destructor deallocates memory allocated by the contsructor
!     -----------------------------------------------------------------
!
      SUBROUTINE DestructKeyValueArray(this)
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(KeyValueArray) :: this
         
      DEALLOCATE (this%items)
      this%constructed = .FALSE.
      this%size        = 0
      this%noOfItems   = 0
!
      END SUBROUTINE DestructKeyValueArray
!
!     /////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------------
!!    KeyValueArrayIsConstructed returns true if the constructor
!!    has been called for "this" instance. This is useful for avoiding
!!    errors but is not always necessary.
!     -------------------------------------------------------------------------
!
      LOGICAL FUNCTION KeyValueArrayIsConstructed(this)
!
        TYPE(KeyValueArray) :: this
        KeyValueArrayIsConstructed = this%constructed
      
      END FUNCTION KeyValueArrayIsConstructed
!
!     /////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------------
!!    This routine is used to add an item to the list. If the list is
!!    too small, then it is reallocated with an additional KV_ALLOC_SIZE
!!    number of items.
!     -------------------------------------------------------------------------
!
      SUBROUTINE AddItemToKeyValueArray( item, this )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(KeyValueArray) :: this
      TYPE(KeyValuePair) :: item
!      
      IF ( this%noOfItems == 0 )     THEN
         this%items( 1 ) = item
         this%noOfItems = 1
      ELSE
         this%noOfItems = this%noOfItems + 1
         IF ( this%noOfItems > this%size ) CALL ReSize( this )
         this%items( this%noOfItems )   = item
      END IF
!      
      END SUBROUTINE AddItemToKeyValueArray
!
!     /////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------
!!    Returns the key-value pair that has the key value.
!     -----------------------------------------------------------------
!
      TYPE(KeyValuePair) FUNCTION ItemForKeyInList( key, this ) RESULT(pair)
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(KeyValueArray) :: this
      CHARACTER( LEN=* ) :: key
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER            :: j
!
      CALL SetPair( pair, "undefined", "undefined" )
      
      DO j = 1, this%noOfItems ! mostly should be the first one
      
         IF ( TRIM(this%items(j)%key) == TRIM(key) )     THEN
            pair = this%items(j)
            RETURN
         END IF
         
      END DO
      
      END FUNCTION ItemForKeyInList
!
!     /////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------
!!    ReSize reallocates memory for the keyValue List. This is private.
!!    Deconstructs the list if the resize fails.
!     -----------------------------------------------------------------
!
      SUBROUTINE ReSize( this )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE( KeyValueArray ) :: this
!
!     ---------------
!     Local Variables
!     ---------------
!
      TYPE(KeyValuePair), DIMENSION( : ), ALLOCATABLE :: tempList
      INTEGER                                         :: k, oldSize, status
!
!     -----------------------------
!     Make temparary copy of arrays
!     -----------------------------
!
      ALLOCATE( tempList( this%size ), STAT=status )
      IF ( status /= 0 )     THEN
         PRINT *,  "KeyValueArray ReSize: Allocation failure"
         RETURN
      END IF
      DO k = 1, this%size
         tempList( k ) = this%items( k )
      END DO
!
!     --------------------
!     Delete current items
!     --------------------
!
      DEALLOCATE( this%items )
!
!     ----------------------
!     Allocate and copy back
!     ----------------------
!
      oldSize    = this%size
      this%size  = this%size + KV_ALLOC_SIZE
      ALLOCATE( this%items( this%size ), STAT=status )
      IF ( status /= 0 )     THEN
         PRINT *,  "KeyValueArray ReSize: Allocation failure"
      ELSE
         DO k = 1, oldSize
            this%items( k ) = tempList( k )
         END DO
      END IF
!
!     -------
!     Cleanup
!     -------
!
      DEALLOCATE( tempList )
!
      END SUBROUTINE ReSize
!
!
!                  ************************************
!                  Iterator Methods for KeyValueArray's
!                  ************************************
!
!     /////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    Iterator Method: InitDictionaryIterator sets the dictionary to
!!    the first entry in the dictionary.
!     ----------------------------------------------------------------
!
      SUBROUTINE InitKeyValueArrayIterator(this)
         TYPE( KeyValueArray ) :: this
         this%iteratorItemID = 0
      END SUBROUTINE InitKeyValueArrayIterator
!
!     /////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------
!!    Iterator Method: InitKeyValueArrayIteratorNext Returns the next item in 
!!    the dictionary.
!     -----------------------------------------------------------------
!
      FUNCTION KeyValueArrayIteratorNext(this) RESULT(p)
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(KeyValueArray) :: this
!
!     ---------------
!     Local Variables
!     ---------------
!
      TYPE(KeyValuePair), POINTER :: p
      INTEGER                     :: start, end, j
!
!     --------------------------------------------------
!     Start at current marker and run to end of the list
!     --------------------------------------------------
!
      start = this%iteratorItemID + 1
      end   = this%size
      DO j = start, end
         this%iteratorItemID = j
         p => this%items(j)
         IF ( .NOT.Undefined(p) )     RETURN
      END DO
!
!     -------------------------------
!     End of list hit, return nothing 
!     -------------------------------
!
      NULLIFY (p)

      END FUNCTION KeyValueArrayIteratorNext
!
!     /////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------
!!    PrintKeyValueArray is a debugging routine that prints the key/value
!!    pairs to the specified file unit
!     -----------------------------------------------------------------
!
      SUBROUTINE PrintKeyValueArray( this, unit )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE( KeyValueArray )  :: this
      INTEGER               :: unit
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER :: j
      
      DO j = 1, this%noOfItems
         IF ( .NOT.Undefined(this%items(j)) )   THEN
            WRITE( unit, * ) TRIM(this%items(j)%key), " = ", &
     &                       TRIM(this%items(j)%value)
         END IF
      END DO
      
      END SUBROUTINE PrintKeyValueArray
!
!  **********     
   END MODULE KeyValueArrayClass
!  **********     
!
! /////////////////////////////////////////////////////////////////////
!
!
!     DictionaryClass
!
!!     Created on: June 14, 2002
!!
!!     Modification History:
!!        June 18, 2002 : Add reSize and error handler
!!        Nov 26, 2002  : Use a simple Hash table for lookup
!!
!!     This class defines a mutable dictionary. To the user it defines
!      a keyword/value relationship. A dictionary has an error handler
!!     instance variable that can be polled to test the success of a
!!     call. A dictionary can be one of two types: CASE_SENSITIVE_DICTIONARY
!!     and CASE_INSENSITIVE_DICTIONARY depending on of the user wants
!!     to care if the key is case sensitive. The default is 
!!     CASE_SENSITIVE_DICTIONARY.
!
!      TYPE Dictionary
!
!
!      PUBLIC METHODS:
!
!         Constructor/Destructor:
!
!            SUBROUTINE ConstructDictionary( this )
!            SUBROUTINE ConstructDictionary( this, type )
!            SUBROUTINE DestructDictionary ( this )
!
!         Adding entries to dictionary:
!
!            SUBROUTINE AddValue_ForKey_ToDict_( value, key, this )
!
!         Getting values from dictionary:
!
!            SUBROUTINE GetValue_ForKey_FromDict_( val, key, dictionary )
!               val can be INTEGER, REAL( RP ), CHARACTER or LOGICAL
!
!         Debugging
!
!            SUBROUTINE PrintDictionary( this, unit )
!         
!         Inquiry:
!
!            LOGICAL FUNCTION DictionaryIsConstructed( this )
!            LOGICAL FUNCTION KeyFoundInDictionary( key, dictionary )
!            INTEGER FUNCTION NumberOfEntriesInDictionary( dictionary )
!            TYPE(KeyValuePair) FUNCTION KeyValuePairForKeyFromDict( key, this )
!
!         Iterator Methods:
!
!            SUBROUTINE InitDictionaryIterator(this)
!            TYPE(KeyValuePair), POINTER FUNCTION DictionaryIteratorNext(this)
!
!      PRIVATE METHODS:
!
!         INTEGER FUNCTION Hash(s,n)
!                     
!!
!!     Usage:
!!        To create a dictionary, call the constructor. When it is
!!        no longer needed, call the destructor. To add keyword pairs, 
!!        call the AddValue method. Note that, currently, the value
!!        is assumed to stored as a string of length DICT_VALUE_STRING_LENGTH, 
!!        so if computed values are to be stored in a dictionary, they
!!        need to be converted to a string first before being added. Finally, 
!!        to access a dictionary entry, use GetValue_ForKey_FromDict_. One can
!!        Iterate through a dictionary to find all keys/values using the
!!        iterator methods.
!!
!!     EXAMPLE:
!!
!!>    TYPE( Dictionary )                        :: theDict
!!     CHARACTER( LEN=DICT_KWD_STRING_LENGTH )   :: kw
!!     DOUBLE PRECISION                          :: x
!!     TYPE(KeyValuePair), POINTER               :: p
!!
!!     CALL ConstructDictionary( theDict )
!!     
!!     kw = "aKeyword"
!!     x  = 3.14
!!     CALL AddValue_ForKey_ToDict_( x, kw, theDict )
!!
!!     CALL GetValue_ForKey_FromDict_( x, kw, theDict )
!!     WRITE( 6, * ) x
!!
!!     CALL InitDictionaryIterator(theDict)
!!     DO
!!        p => DictionaryIteratorNext(theDict)
!!        IF ( .NOT.ASSOCIATED(p) )     EXIT
!!        -- DO Something with p --
!!     END DO
!!<
!
!
!!    @author David A. Kopriva
!
! /////////////////////////////////////////////////////////////////////
!
!  ******
   MODULE DictionaryClass
!  ******
!
     USE SMConstants
     USE KeyValueArrayClass
     
     IMPLICIT NONE
     
     PUBLIC  :: Dictionary, GetValue_ForKey_FromDict_, NumberOfEntriesInDictionary
     PRIVATE :: ReSize
!
!    ----------
!    Interface:
!    ----------
!
!    -----------------------------------------
!!   Public interface to the dictionary values
!    -----------------------------------------
!
     INTERFACE GetValue_ForKey_FromDict_
        MODULE PROCEDURE GetIntForKeyFromDict
        MODULE PROCEDURE GetLogicalForKeyFromDict
        MODULE PROCEDURE GetStringForKeyFromDict
        MODULE PROCEDURE GetRealForKeyFromDict
     END INTERFACE
!
!    -----------------------------------
!!   Public interface to the Add methods
!    -----------------------------------
!
     INTERFACE AddValue_ForKey_ToDict_
        MODULE PROCEDURE AddIntValueForKeyToDict
        MODULE PROCEDURE AddRealValueForKeyToDict
        MODULE PROCEDURE AddLogicalValueForKeyToDict
        MODULE PROCEDURE AddStringValueForKeyToDict
     END INTERFACE
!
!    -----------------
!    Module constants:
!    -----------------
!!
!       -----------------------------------------------------------
!!      The allocation size determines the number of memory chunks
!!      when memory is allocated.
!       -----------------------------------------------------------
!
        INTEGER, PARAMETER, PUBLIC  :: CASE_SENSITIVE_DICTIONARY   = 0
        INTEGER, PARAMETER, PUBLIC  :: CASE_INSENSITIVE_DICTIONARY = 1
        INTEGER, PARAMETER, PRIVATE :: DICT_ALLOC_SIZE             = 51
!
!    ----------------
!    Type definitions
!    ----------------
!
!       ----------------------------------------------------------------
!!      A dictionary defines keyword based access to data. In this case
!!      there is a correspondence between a keyword and some typed data.
!!      The typed data is stored as a string and then converted as
!!      necessary to allow for multiple data types to be stored. However, 
!!      this fact need not concern the user.
!!
!       ----------------------------------------------------------------
!
     TYPE Dictionary
        LOGICAL                                      :: constructed
        INTEGER                                      :: size
        INTEGER                                      :: numEntries
        INTEGER                                      :: dictionaryType
        TYPE( KeyValueArray ), DIMENSION(:), POINTER :: entries
        INTEGER                                      :: iteratorEntryID
     END TYPE Dictionary
!
!    ========
     CONTAINS
!    ========
!
!     /////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    Constructor: Allocate memory for the dictionary
!     ----------------------------------------------------------------
!
      SUBROUTINE ConstructDictionary( this, type )
!
      TYPE( Dictionary )                  :: this
      INTEGER, OPTIONAL                   :: type
      
      INTEGER                             :: status
!
      INTEGER :: j
!
      IF ( PRESENT(type) )     THEN
         this%dictionaryType = type
      ELSE
         this%dictionaryType = CASE_SENSITIVE_DICTIONARY
      END IF
            
      this%numEntries = 0
      this%size       = DICT_ALLOC_SIZE
!
!     ----------------
!     Allocate Storage
!     ----------------
!      
      ALLOCATE( this%entries( DICT_ALLOC_SIZE ), STAT=status )
               
      IF ( status == 0 )     THEN
         this%constructed = .TRUE.
      ELSE
         PRINT *, "ConstructDictionary: Allocation failure"
         RETURN
      END IF
!
!     ------------------------------------------
!     Construct the entries with nothing in them
!     ------------------------------------------
!
      DO j = 1, this%size
         CALL ConstructKeyValueArray(this%entries(j))
      END DO
!
      END SUBROUTINE ConstructDictionary
!
!     /////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    Destructor:: Clears allocated memory
!     ----------------------------------------------------------------
!
      SUBROUTINE DestructDictionary( this )
!
      TYPE( Dictionary ) :: this
      INTEGER            :: j
      
      DO j = 1, this%size
         CALL DestructKeyValueArray( this%entries(j) )
      END DO

      this%numEntries = 0
      this%size       = 0
      this%constructed = .FALSE.
      
      END SUBROUTINE DestructDictionary
!
!     /////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    IsConstructed returns true if the constructor has been called
!!    for "this"
!     ----------------------------------------------------------------
!
      LOGICAL FUNCTION DictionaryIsConstructed( this )
!
         TYPE( Dictionary ) :: this
         DictionaryIsConstructed = this%constructed
      
      END FUNCTION DictionaryIsConstructed
!
!     /////////////////////////////////////////////////////////////////
!
!     ------------------------------------------------------
!!    KeyFoundInDictionary returns .TRUE. if the key exists
!!    in the dictionary.
!     ------------------------------------------------------
!
      LOGICAL FUNCTION KeyFoundInDictionary( key, this )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE( Dictionary ) :: this
      CHARACTER( LEN=* ) :: key
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER                 :: j
      TYPE(KeyValuePair)      :: pair
      CHARACTER(LEN=LEN(key)) :: s
      
      IF ( this%dictionaryType == CASE_INSENSITIVE_DICTIONARY )     THEN
         j    = Hash( s, this%size )
         pair = ItemForKeyInList( s, this%entries(j) )
      ELSE
         j    = Hash( key, this%size )
         pair = ItemForKeyInList( key, this%entries(j) )
      END IF
      
      IF ( Undefined(pair) )     THEN
         KeyFoundInDictionary = .FALSE.
      ELSE
         KeyFoundInDictionary = .TRUE.
      END IF
     
     END FUNCTION KeyFoundInDictionary
!
!     /////////////////////////////////////////////////////////////////
!
!     ------------------------------------------------------
!!    KeyFoundInDictionary returns .TRUE. if the key exists
!!    in the dictionary.
!     ------------------------------------------------------
!
      INTEGER FUNCTION NumberOfEntriesInDictionary( this )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE( Dictionary ) :: this
      
      NumberOfEntriesInDictionary = this%numEntries
     
     END FUNCTION NumberOfEntriesInDictionary
!
!     /////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------
!!    ReSize reallocates memory for the dictionary. This is private.
!     -----------------------------------------------------------------
!
      SUBROUTINE ReSize( this )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE( Dictionary ) :: this
!
!     ---------------
!     Local Variables
!     ---------------
!
      TYPE(KeyValueArray), DIMENSION(:), ALLOCATABLE :: tempList
      INTEGER                                       :: k, oldSize, status
!
!     -----------------------------
!     Make temparary copy of arrays
!     -----------------------------
!
      ALLOCATE( tempList( this%size ), STAT=status )
      IF ( status /= 0 )     THEN
         PRINT *,  "Dictionary ReSize: Allocation failure"
         RETURN
      END IF
      DO k = 1, this%size
         tempList( k ) = this%entries( k )
      END DO
!
!     ---------------------
!     Delete current arrays
!     ---------------------
!
      DEALLOCATE( this%entries )
!
!     ----------------------
!     Allocate and copy back
!     ----------------------
!
      oldSize    = this%size
      this%size  = this%size + DICT_ALLOC_SIZE
      ALLOCATE( this%entries( this%size ), STAT=status )
      IF ( status /= 0 )     THEN
         PRINT *,  "Dictionary ReSize: Allocation failure"
      ELSE
         DO k = 1, oldSize
            this%entries( k ) = tempList( k )
         END DO
      END IF
!
!     -------
!     Cleanup
!     -------
!
      DEALLOCATE( tempList )
!
      END SUBROUTINE ReSize
!
!     //////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------------
!!    AddStringValueForKeyToDict inserts the key/value pair into the dictionary
!     -------------------------------------------------------------------------
!
      SUBROUTINE AddStringValueForKeyToDict( value, key, this )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE( Dictionary ) :: this
      CHARACTER( LEN=* ) :: key
      CHARACTER( LEN=* ) :: value
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER                 :: j
      TYPE(KeyValuePair)      :: pair
      CHARACTER(LEN=LEN(key)) :: s
      
      IF ( this%dictionaryType == CASE_INSENSITIVE_DICTIONARY )     THEN
         j    = Hash( s, this%size )
         CALL SetPair( pair, s, value )
      ELSE
         j    = Hash( key, this%size )
         CALL SetPair( pair, key, value )
      END IF

      CALL AddItemToKeyValueArray( pair, this%entries(j) )
      this%numEntries = this%numEntries + 1
      
      END SUBROUTINE AddStringValueForKeyToDict
!
!     //////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------------
!!    AddIntValueForKeyToDict inserts the key/int value pair into the dictionary
!     -------------------------------------------------------------------------
!
      SUBROUTINE AddIntValueForKeyToDict( value, key, this )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE( Dictionary ) :: this
      INTEGER            :: value
      CHARACTER( LEN=* ) :: key
!
!     ---------------
!     Local Variables
!     ---------------
!
      CHARACTER( LEN=DICT_VALUE_STRING_LENGTH ) :: strValue
      
      WRITE(strValue,*) value
      CALL AddStringValueForKeyToDict( strValue, key, this )
     
      END SUBROUTINE AddIntValueForKeyToDict
!
!     //////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------------
!!    AddLogicalValueForKeyToDict inserts the key/int value pair into the 
!!    dictionary
!     -------------------------------------------------------------------------
!
      SUBROUTINE AddLogicalValueForKeyToDict( value, key, this )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE( Dictionary ) :: this
      LOGICAL            :: value
      CHARACTER( LEN=* ) :: key
!
!     ---------------
!     Local Variables
!     ---------------
!
      CHARACTER( LEN=DICT_VALUE_STRING_LENGTH ) :: strValue
      
      WRITE(strValue,*) value
      CALL AddStringValueForKeyToDict( strValue, key, this )
     
      END SUBROUTINE AddLogicalValueForKeyToDict
!
!     //////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------------
!!    AddRealValueForKeyToDict inserts the key/int value pair into the 
!!    dictionary
!     -------------------------------------------------------------------------
!
      SUBROUTINE AddRealValueForKeyToDict( value, key, this )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE( Dictionary ) :: this
      REAL(KIND=RP)  :: value
      CHARACTER( LEN=* ) :: key
!
!     ---------------
!     Local Variables
!     ---------------
!
      CHARACTER( LEN=DICT_VALUE_STRING_LENGTH ) :: strValue
      
      WRITE(strValue,*) value
      CALL AddStringValueForKeyToDict( strValue, key, this )
     
      END SUBROUTINE AddRealValueForKeyToDict
!
!     //////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------------
!!    KeyValuePairForKeyFromDict finds the key/value pair for a particular key
!     -------------------------------------------------------------------------
!
      TYPE(KeyValuePair) FUNCTION KeyValuePairForKeyFromDict( key, this ) &
     &                                                        RESULT(pair)
!
!     ---------
!     Arguments
!     ---------
!
      TYPE( Dictionary ) :: this
      CHARACTER( LEN=* ) :: key
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER                 :: j
      CHARACTER(LEN=LEN(key)) :: s
      
      IF ( this%dictionaryType == CASE_INSENSITIVE_DICTIONARY )     THEN
         j    = Hash( s, this%size )
         pair = ItemForKeyInList( s, this%entries(j) )
      ELSE
         j    = Hash( key, this%size )
         pair = ItemForKeyInList( key, this%entries(j) )
      END IF
      
      END FUNCTION KeyValuePairForKeyFromDict
!
!     /////////////////////////////////////////////////////////////////
!
!     ------------------------------------------------------
!!    GetIntForKeyFromDict returns an integer value associated 
!!    with the given key.
!     ------------------------------------------------------
!
      SUBROUTINE GetIntForKeyFromDict( val, key, this )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE( Dictionary ) :: this
      CHARACTER( LEN=* ) :: key
      INTEGER          :: val
!
!     ---------------
!     Local Variables
!     ---------------
!
      CHARACTER( LEN=LINE_LENGTH ) :: message
      TYPE(KeyValuePair)                 :: pair
     
      pair = KeyValuePairForKeyFromDict( key, this )
      
      IF ( Undefined(pair) )     THEN
        message = "Keyword '" // TRIM( key ) // "' not found"
        PRINT *, message
      ELSE
         IF( INDEX(pair%value,".") /= 0 )     THEN
            message = "Non-integer value found for Keyword " // TRIM( key )
            PRINT *, message
         ELSE
           READ( pair%value, * ) val
         END IF
      END IF
     
     END SUBROUTINE GetIntForKeyFromDict
!
!     /////////////////////////////////////////////////////////////////
!
!     ------------------------------------------------------
!!    GetLogicalForKey returns a logical value associated
!!    with the given key.
!     ------------------------------------------------------
!
      SUBROUTINE GetLogicalForKeyFromDict( val, key, this )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE( Dictionary ) :: this
      CHARACTER( LEN=* ) :: key
      LOGICAL          :: val
!
!     ---------------
!     Local Variables
!     ---------------
!
      CHARACTER( LEN=LINE_LENGTH ) :: message
      TYPE(KeyValuePair)                 :: pair
     
      pair = KeyValuePairForKeyFromDict( key, this )
      
      IF ( Undefined(pair) )     THEN
        message = "Keyword '" // TRIM( key ) // "' not found"
        PRINT *, message
      ELSE
        READ( pair%value, * ) val
      END IF
     
      END SUBROUTINE GetLogicalForKeyFromDict
!
!     /////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------
!!    GetStringForKeyFromDict returns an string value associated 
!!    with the given key.
!     ----------------------------------------------------------
!
      SUBROUTINE GetStringForKeyFromDict( val, key, this )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE( Dictionary ) :: this
      CHARACTER( LEN=* ) :: key
      CHARACTER( LEN=* ) :: val
!
!     ---------------
!     Local Variables
!     ---------------
!
      CHARACTER( LEN=LINE_LENGTH ) :: message
      TYPE(KeyValuePair)                 :: pair
     
      pair = KeyValuePairForKeyFromDict( key, this )
      
      IF ( Undefined(pair) )     THEN
        message = "Keyword '" // TRIM( key ) // "' not found"
        PRINT *, message
        val = ""
      ELSE
        val = pair%value
      END IF
     
      END SUBROUTINE GetStringForKeyFromDict
!
!     /////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------
!!    GetRealForKeyFromDict returns a real kind=RP value associated 
!!    with the given key.
!     -----------------------------------------------------------------
!
      SUBROUTINE GetRealForKeyFromDict( val, key, this )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE( Dictionary )  :: this
      CHARACTER( LEN=* )  :: key
      REAL( KIND=RP ) :: val
!
!     ---------------
!     Local Variables
!     ---------------
!
      CHARACTER( LEN=LINE_LENGTH ) :: message
      TYPE(KeyValuePair)                 :: pair
     
      pair = KeyValuePairForKeyFromDict( key, this )
      
      IF ( Undefined(pair) )     THEN
        message = "Keyword '" // TRIM( key ) // "' not found"
        PRINT *, message
      ELSE
        READ( pair%value, * ) val
      END IF
     
      END SUBROUTINE GetRealForKeyFromDict
!
! /////////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------
!!    This routine returns the Hash value for a string. This is not the
!!    most sophisticated function I've found, but it is simple.
!     -----------------------------------------------------------------
!
      INTEGER FUNCTION Hash(s,n)
      
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      CHARACTER(LEN=*):: s !! The string to Hash
      INTEGER         :: n !! The size of the Hash table
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER         :: h !! Intermediate value of Hash
      INTEGER         :: j
      
      h = 0
      DO j = 1,LEN_TRIM(s)
         h = h + ICHAR(s(j:j))
      END DO
      Hash = MOD(h,n) + 1
      
      END FUNCTION Hash
!
!     /////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    Iterator Method: InitDictionaryIterator sets the dictionary to
!!    the first entry in the dictionary.
!     ----------------------------------------------------------------
!
      SUBROUTINE InitDictionaryIterator(this)
         TYPE( Dictionary ) :: this
         this%iteratorEntryID = 1
         CALL InitKeyValueArrayIterator( this%entries(this%iteratorEntryID) )
      END SUBROUTINE InitDictionaryIterator
!
!     /////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------
!!    Iterator Method: DictionaryIteratorNext Returns the next item in 
!!    the dictionary.
!     -----------------------------------------------------------------
!
      FUNCTION DictionaryIteratorNext(this) RESULT(p)
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(Dictionary) :: this
!
!     ---------------
!     Local Variables
!     ---------------
!
      TYPE(KeyValuePair), POINTER :: p
      INTEGER                     :: start, end, j
!
      p => KeyValueArrayIteratorNext( this%entries(this%iteratorEntryID) )
      IF ( ASSOCIATED(p) )     RETURN
!
!     ------------------------------------------------------------------
!     Nothing (more) found in this KeyValueArray, so continue through the
!     entries until one is found
!     ------------------------------------------------------------------
!
      start = this%iteratorEntryID + 1
      end = this%size
      DO j = start, end
         this%iteratorEntryID = j
         CALL InitKeyValueArrayIterator( this%entries(j) )
         p => KeyValueArrayIteratorNext( this%entries(j) )
         IF ( ASSOCIATED(p) )     RETURN
      END DO
      NULLIFY (p)
      
      END FUNCTION DictionaryIteratorNext
!
!     /////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------
!!    PrintDictionary is a debugging routine that prints the key/value
!!    pairs to the specified file unit
!     -----------------------------------------------------------------
!
      SUBROUTINE PrintDictionary( this, unit )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE( Dictionary )  :: this
      INTEGER             :: unit
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER :: j
      
      WRITE( unit, * ) "Key/Value Pairs in dictionary:"
      WRITE( unit, * ) " "
      DO j = 1, this%size
         CALL PrintKeyValueArray( this%entries(j), unit )
      END DO
      
      END SUBROUTINE PrintDictionary
!
!  **********     
   END MODULE DictionaryClass
!  **********     
