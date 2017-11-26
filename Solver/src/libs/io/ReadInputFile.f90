!
!////////////////////////////////////////////////////////////////////////
!
!      ReadInputFile.f90
!      Created: June 10, 2015 at 3:09 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
! 
      SUBROUTINE ReadInputFile (controlVariables)
         USE SMConstants
         USE FTValueDictionaryClass
         USE SharedBCModule
         USE mainKeywordsModule
         use MPI_Process_Info 
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(FTValueDictionary) :: controlVariables
!
!        ---------------
!        Local variables
!        ---------------
!
         CHARACTER(LEN=LINE_LENGTH) :: inputLine
         CHARACTER(LEN=LINE_LENGTH) :: keyword, keywordValue
         CHARACTER(LEN=LINE_LENGTH) :: boundaryName
         CHARACTER(LEN=LINE_LENGTH) :: boundaryType
         CHARACTER(LEN=LINE_LENGTH) :: boundaryValue
         CHARACTER(LEN=LINE_LENGTH) :: arg
         character(len=LINE_LENGTH) :: boundaryNameControlVariable
         INTEGER                    :: numberOfBCs, k
         INTEGER                    :: ist
         logical                                 :: isInsideHagstagZone
!
!        ---------------------------------------
!        External functions from FileReading.f90
!        ---------------------------------------
!
         REAL(KIND=RP)             , EXTERNAL    :: GetRealValue
         INTEGER                   , EXTERNAL    :: GetIntValue
         CHARACTER(LEN=LINE_LENGTH), EXTERNAL    :: GetStringValue, GetKeyword, GetValueAsString
         LOGICAL                   , EXTERNAL    :: GetLogicalValue
         
!
!        -----------------------------------------------
!        Read the input file.
!
!        we use dictionaries to store the input file 
!        parameters.
!        -----------------------------------------------
!

         if ( command_argument_count() .eq. 0 ) then
            if ( MPI_Process % isRoot ) then
               write(STD_OUT,'(/,/,A)') "*** ERROR: Missing input file"
               write(STD_OUT,'(A)') "*** ERROR: Syntax is: "
               write(STD_OUT,'(A)') "*** ERROR:             >> HORSES3D ControlFile.control"
               write(STD_OUT,'(A,/,/)') "*** ERROR: Stopping."
            end if
            stop
         end if

         CALL get_command_argument(1, arg)
         OPEN(UNIT=10,FILE=trim(arg))

         isInsideHagstagZone = .false.

         DO
            READ(10,'(A132)', IOSTAT = ist) inputLine
            IF(ist /= 0 ) EXIT !.OR. inputLine(1:1) == '/'
            IF ( inputLine(1:1) == "!" .OR. inputLine(1:1) == '/') CYCLE ! Skip comments

            if ( index(inputLine,'#define') .ne. 0 ) then
               isInsideHagstagZone = .true.

            elseif ( (index(inputLine,'#end') .ne. 0) .and. (isInsideHagstagZone) ) then
               isInsideHagstagZone = .false.

            end if

            if ( isInsideHagstagZone ) cycle
            
            keyword      = ADJUSTL(GetKeyword(inputLine))
            keywordValue = ADJUSTL(GetValueAsString(inputLine))
            CALL toLower(keyword)
            CALL controlVariables % addValueForKey(keywordValue,TRIM(keyword))
            
            IF(keyword == numberOfBoundariesKey) THEN 
!
!              ---------------------------------------------------------------------------
!              We will store the type and values of the boundaries in dictionaries so that
!              we can associate a name of a boundary curve found in the mesh file with a
!              particular value and type of boundary conditions.
!              ---------------------------------------------------------------------------
!
               numberOfBCs = controlVariables%integerValueForKey(numberOfBoundariesKey)
               
               DO k = 1, numberOfBCs 
                  READ(10,*) boundaryName, boundaryValue, boundaryType
                  CALL toLower(boundaryName)
                  CALL toLower(boundaryType)
                  CALL bcTypeDictionary % addValueForKey(boundaryType, boundaryName)
                  CALL bcValueDictionary % addValueForKey(boundaryValue, boundaryName)
                  write(boundaryNameControlVariable,'(A,I0)') "BoundaryName",k
                  call controlVariables % addValueForKey(boundaryName,trim(boundaryNameControlVariable))
               END DO
            END IF
            
            IF(keyword == "boundaries to analyze") THEN 
!
!              ---------------------------------------------------------------------------
!              We will store the name of the zones (boundaries where we want to compute lift and
!              drag) in a dictionary. Currently, this quantities are computed taking into 
!              account the angle of attack. However, this dictionary can be easily modified in 
!              order to compute forces in a different direction (added as value).
!              ---------------------------------------------------------------------------
!
               numberOfBCs = controlVariables%integerValueForKey("boundaries to analyze")
               
               DO k = 1, numberOfBCs 
                  READ(10,*) boundaryName
                  
                  CALL toLower(boundaryName)
                  CALL zoneNameDictionary % addValueForKey(boundaryName, boundaryName)
                  
               END DO
            END IF
            
         END DO

         CLOSE(UNIT=10)

      END SUBROUTINE ReadInputFile
      
