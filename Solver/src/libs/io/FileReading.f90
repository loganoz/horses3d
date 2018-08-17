!
!////////////////////////////////////////////////////////////////////////
!
!      FileReading.f90
!      Created: 2010-08-27 12:14:31 -0400 
!      By: David Kopriva  
!
!///////////////////////////////////////////////////////////////////////
!
module FileReadingUtilities
   USE SMConstants
   use RealDataLinkedList
   implicit none
   
   public
   
contains
!
!///////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    Extracts the string at the left of the =, which corresponds to
!!    the keyword for an input value.
!     ----------------------------------------------------------------
!
      CHARACTER( LEN=LINE_LENGTH ) FUNCTION GetKeyword( inputLine )
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
!!    Extracts the string within the quotes in an input file
!     ----------------------------------------------------------------
!
      CHARACTER( LEN=LINE_LENGTH ) FUNCTION GetValueAsString( inputLine )
         IMPLICIT NONE
!
         CHARACTER ( LEN = * ) :: inputLine
         INTEGER               :: strLen, leq
         INTEGER               :: cStart, cEnd
!
         cStart = INDEX(inputLine,'"')
         IF ( cStart /= 0 )     THEN
             cEnd             = INDEX(inputLine, '"', .true. )
             GetValueAsString = inputLine( cStart+1: cEnd-1 )
         ELSE 
            leq              = INDEX( inputLine, '=' )
            strLen           = LEN_TRIM( inputLine )
            GetValueAsString = inputLine( leq+1: strLen )
         END IF 
!
!
      END FUNCTION GetValueAsString
!
!///////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    "Read" the "value" of an real number declared
!!     after an = sign in an inputLine
!     ----------------------------------------------------------------
!
      INTEGER FUNCTION GetRealValue( inputLine )
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

      character(len=LINE_LENGTH) function RemovePath( inputLine )
         implicit none
         character(len=*)     :: inputLine
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: pos
!
!        Get the last forward slash ocurrence
!        ------------------------------------
         pos = index(inputLine,'/',BACK=.true.)

         if ( pos .eq. 0 ) then
            RemovePath = inputLine

         else
            RemovePath = inputLine(pos+1:)

         end if
         
      end function RemovePath

      character(len=LINE_LENGTH) function getFileName( inputLine )
         implicit none
         character(len=*)     :: inputLine
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: pos
!
!        Get the last point ocurrence
!        ----------------------------
         pos = index(inputLine,'.',BACK=.true.)

         if ( pos .eq. 0 ) then  
            getFileName = inputLine

         else
            getFileName = inputLine(1:pos-1)

         end if

      end function getFileName
      
      character(len=LINE_LENGTH) function getFileExtension( inputLine )
         implicit none
         character(len=*)     :: inputLine
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: pos
!
!        Get the last point ocurrence
!        ----------------------------
         pos = index(inputLine,'.',BACK=.true.)

         if ( pos .eq. 0 ) then  
            getFileExtension = inputLine

         else
            getFileExtension = inputLine(pos+1:)

         end if

      end function getFileExtension

      function getRealArrayFromString( line ) result ( array )
!
!           ****************************************************
!                    Gets an array from a string of the 
!              form: 
!                       line = "[a,b,c,...]"
!           ****************************************************
!
         implicit none
         character(len=*),    intent(in)  :: line
         real(kind=RP), allocatable       :: array(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: pos1 , pos2 , pos
         character(len=LINE_LENGTH)   :: auxline
         type(RealDataLinkedList_t)    :: Data
         real(kind=RP)                 :: value
         integer  :: io
        
         pos1 = index(line,"[")
         pos2 = index(line,"]") 

         if ( (pos1 .eq. 0) .or. (pos2 .eq. 0) ) then
!
!           There are no brackets in the string
!           -----------------------------------
            return
         end if
         
         auxline = line(pos1+1:pos2-1)
!
!        Get the elements
!        ----------------
         do
            pos = index(auxline , "," ) 

            if ( pos .gt. 0 ) then

               read(auxline(1:pos-1),*,iostat=io) value 
               if ( io .lt. 0 ) then
                  return
               end if

               call Data % Append(value)
   
               auxline = auxline(pos+1:)

            else
               read(auxline ,*,iostat=io) value 
               if ( io .lt. 0 ) then
                  return
               end if

               call Data % append(value)

               exit

            end if

         end do

         call Data % load(array)
         call Data % destruct

      end function getRealArrayFromString
      
      function getCharArrayFromString( line , linelength) result ( array )
!
!           ****************************************************
!                    Gets a character array from a string of the 
!              form: 
!                       line = "[a,b,c,...]"
!           ****************************************************
!
         implicit none
         character(len=*)       ,    intent(in)  :: line
         integer                                 :: linelength
         character(len=linelength), allocatable  :: array(:)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: pos1 , pos2 , pos
         character(len=LINE_LENGTH) :: auxline
         integer                    :: io
         integer                    :: numOfElems, elem
         
         numOfElems = 1
         
         pos1 = index(line,"[")
         pos2 = index(line,"]") 

         if ( (pos1 .eq. 0) .or. (pos2 .eq. 0) ) then
!
!           There are no brackets in the string
!           -----------------------------------
            return
         end if
         
!
!        Get number of elements
!        ----------------------
         auxline = line(pos1+1:pos2-1)
         do
            pos = index(auxline , "," ) 

            if ( pos .gt. 0 ) then

               numOfElems = numOfElems + 1
               auxline = auxline(pos+1:)
            else
               exit
            end if

         end do
         
         allocate ( array(numOfElems) )

!
!        Get the elements
!        ----------------
         auxline = line(pos1+1:pos2-1)
         elem = 1
         do
            pos = index(auxline , "," ) 

            if ( pos .gt. 0 ) then

               read(auxline(1:pos-1),*,iostat=io) array(elem) 
               if ( io .lt. 0 ) then
                  return
               end if
   
               auxline = auxline(pos+1:)
               elem    = elem + 1

            else
               read(auxline ,*,iostat=io) array(elem)  
               if ( io .lt. 0 ) then
                  return
               end if

               exit

            end if

         end do

      end function getCharArrayFromString
end module FileReadingUtilities
