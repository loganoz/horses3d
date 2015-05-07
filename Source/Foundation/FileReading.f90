!
!////////////////////////////////////////////////////////////////////////
!
!      FileReading.f90
!      Created: 2010-08-27 12:14:31 -0400 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    "Read" the "value" of a real number declared
!!     after an = sign in a inputLine
!     ----------------------------------------------------------------
!
      FUNCTION GetRealValue( inputLine )
         USE SMConstants
         IMPLICIT NONE 
!
         CHARACTER(len = *) :: inputLine
         REAL( KIND=RP )    :: value, GetRealValue
         INTEGER            :: strLen, leq
!
         leq    = INDEX( inputLine, '=' )
         strLen = LEN_TRIM( inputLine )
         READ( inputLine( leq+1:strLen ), * ) value
         GetRealValue = value
!
      END FUNCTION GetRealValue
!
!///////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    "Read" the "value" of an real array declared
!!     after an = sign in an inputLine, e.g. rArray(1:2) = [3.0,4.0]
!     ----------------------------------------------------------------
!
      FUNCTION GetRealArray( inputLine ) RESULT(x)
         USE SMConstants
         IMPLICIT NONE
!
         REAL(KIND=RP), DIMENSION(2) :: x
         
         CHARACTER ( LEN = * ) :: inputLine
         INTEGER               :: cStart, cEnd
!
         cStart = INDEX(inputLine,'[')
         cEnd   = INDEX(inputLine, ']', .true. )
         READ( inputLine( cStart+1: cEnd-1 ), * ) x(1), x(2)
!
      END FUNCTION GetRealArray
!
!///////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    "Read" the "value" of an integer number declared
!!     after an = sign in an inputLine
!     ----------------------------------------------------------------
!
      INTEGER FUNCTION GetIntValue( inputLine )
         IMPLICIT NONE
!
         CHARACTER ( LEN = * ) :: inputLine
         INTEGER               :: value
         INTEGER               :: strLen, leq
!
         leq    = INDEX( inputLine, '=' )
         strLen = LEN_TRIM( inputLine )
         READ( inputLine( leq+1:strLen ), * ) value
         GetIntValue = VALUE
!
      END FUNCTION GetIntValue
!
!///////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    "Read" the "value" of an integer array declared
!!     after an = sign in an inputLine, e.g. iArray(1:2) = [3,4]
!     ----------------------------------------------------------------
!
      FUNCTION GetIntArray( inputLine ) RESULT(N)
         IMPLICIT NONE
!
         INTEGER, DIMENSION(2) :: N
         
         CHARACTER ( LEN = * ) :: inputLine
         INTEGER               :: cStart, cEnd
!
         cStart = INDEX(inputLine,'[')
         cEnd   = INDEX(inputLine, ']', .true. )
         READ( inputLine( cStart+1: cEnd-1 ), * ) N(1), N(2)
!
      END FUNCTION GetIntArray
!
!///////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    Extracts the string within the quotes in an input file
!     ----------------------------------------------------------------
!
      CHARACTER( LEN=LINE_LENGTH ) FUNCTION GetStringValue( inputLine )
         USE SMConstants
         IMPLICIT NONE
!
         CHARACTER ( LEN = * ) :: inputLine
         INTEGER               :: cStart, cEnd
!
         cStart = INDEX(inputLine,'"')
         cEnd   = INDEX(inputLine, '"', .true. )
         GetStringValue = inputLine( cStart+1: cEnd-1 )
!
      END FUNCTION GetStringValue
!
!///////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    "Read" the "value" of an logical declared
!!     after an = sign in an inputLine
!     ----------------------------------------------------------------
!
      LOGICAL FUNCTION GetLogicalValue( inputLine )
         IMPLICIT NONE
!
         CHARACTER ( LEN = * ) :: inputLine
         LOGICAL               :: value
         INTEGER               :: strLen, leq
!
         leq    = INDEX( inputLine, '=' )
         strLen = LEN_TRIM( inputLine )
         READ( inputLine( leq+1:strLen ), * ) value
         GetLogicalValue = value
!
      END FUNCTION GetLogicalValue
