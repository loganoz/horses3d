!
!////////////////////////////////////////////////////////////////////////
!
!      FileReading.f90
!      Created: 2010-08-27 12:14:31 -0400 
!      By: David Kopriva  
!
!///////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    Extracts the string at the left of the ==, which corresponds to
!!    the keyword for an input value.
!     ----------------------------------------------------------------
!
      CHARACTER( LEN=LINE_LENGTH ) FUNCTION GetKeyword( inputLine )
         USE SMConstants
         IMPLICIT NONE
!
         CHARACTER ( LEN = * )          :: inputLine
         CHARACTER ( LEN = LINE_LENGTH) :: adjustedLine
         INTEGER                        :: cEnd
!
         adjustedLine = ADJUSTL(inputLine)
         cEnd         = INDEX(adjustedLine, '=')
         GetKeyword   = adjustedLine(1: cEnd-1 )
!
      END FUNCTION GetKeyword
!
!///////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    "Read" the "value" of an real number declared
!!     after an = sign in an inputLine
!     ----------------------------------------------------------------
!
      INTEGER FUNCTION GetRealValue( inputLine )
         USE SMConstants
         IMPLICIT NONE
!
         CHARACTER ( LEN = * ) :: inputLine
         REAL(KIND=RP)         :: value
         INTEGER               :: strLen, leq
!
         leq    = INDEX( inputLine, '=' )
         strLen = LEN_TRIM( inputLine )
         READ( inputLine( leq+1:strLen ), * ) value
         GetRealValue = VALUE
!
      END FUNCTION GetRealValue
!
!////////////////////////////////////////////////////////////////////////
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
