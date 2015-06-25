!
!////////////////////////////////////////////////////////////////////////
!
!      ReadInputFile.f90
!      Created: June 10, 2015 at 3:09 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
! 
      Module ControlVariablesModule
         USE SMConstants
         IMPLICIT NONE
         
         TYPE NSLiteControlVariables
            CHARACTER(LEN=LINE_LENGTH) :: inputFileName, plotFileName, restartFileName
            REAL(KIND=RP)              :: cfl
            REAL(KIND=RP)              :: tol
            REAL(KIND=RP)              :: mach, RE
            REAL(KIND=RP)              :: AOATheta, AOAPhi
            INTEGER                    :: polynomialOrder
            INTEGER                    :: plotInterval
            INTEGER                    :: numberOfSteps
            INTEGER                    :: numberOfPlotPoints
            LOGICAL                    :: restart
            LOGICAL                    :: flowIsNavierStokes
         END TYPE NSLiteControlVariables
         
      END MODULE ControlVariablesModule
!
!////////////////////////////////////////////////////////////////////////
! 
      SUBROUTINE ReadInputFile (controlVariables)
         USE SMConstants
         USE ControlVariablesModule
         USE SharedBCModule
         
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(NSLiteControlVariables) :: controlVariables
!
!        ---------------
!        Local variables
!        ---------------
!
         CHARACTER(LEN=LINE_LENGTH) :: inputLine
         CHARACTER(LEN=LINE_LENGTH) :: boundaryName
         CHARACTER(LEN=LINE_LENGTH) :: flowEquationsName
         CHARACTER(LEN=LINE_LENGTH) :: boundaryType
         CHARACTER(LEN=LINE_LENGTH) :: boundaryValue
         INTEGER                    :: numberOfBCs, k
!
!        ---------------------------------------
!        External functions from FileReading.f90
!        ---------------------------------------
!
         REAL(KIND=RP)             , EXTERNAL    :: GetRealValue
         INTEGER                   , EXTERNAL    :: GetIntValue
         CHARACTER(LEN=LINE_LENGTH), EXTERNAL    :: GetStringValue
         LOGICAL                   , EXTERNAL    :: GetLogicalValue
!
!        -----------------------------------------------
!        Read the input file.
!
!        Nothing fancy here. Pretty much a fixed format,
!        except that for convenience, we use 
!        dictionaries to store the input file parameters.
!        -----------------------------------------------
!
         READ(5,'(A132)') inputLine
         flowEquationsName = GetStringValue( inputLine )
         
         IF ( flowEquationsName == 'Euler' .OR. flowEquationsName == 'euler' )     THEN
            controlVariables % flowIsNavierStokes = .false.
         ELSE
            controlVariables % flowIsNavierStokes = .true.
         END IF
   
         READ(5,'(A132)') inputLine
         controlVariables % inputFileName = GetStringValue( inputLine )
         
         READ(5,'(A132)') inputLine
         controlVariables % plotFileName = GetStringValue( inputLine )
         
         READ(5,'(A132)') inputLine
         controlVariables % restartFileName = GetStringValue( inputLine )
   
         READ(5,'(A132)') inputLine
         controlVariables % restart = GetLogicalValue( inputLine )
         
         READ(5,'(A132)') inputLine
         controlVariables % polynomialOrder = GetIntValue( inputLine )
         
         READ(5,'(A132)') inputLine
         controlVariables % numberOfSteps = GetIntValue( inputLine )
         
         READ(5,'(A132)') inputLine
         controlVariables % plotInterval = GetIntValue( inputLine )
         
         READ(5,'(A132)') inputLine
         controlVariables % numberOfPlotPoints = GetIntValue( inputLine )
         
         READ(5,'(A132)') inputLine
         controlVariables % tol = GetRealValue( inputLine )
         
         READ(5,'(A132)') inputLine
         controlVariables % cfl = GetRealValue( inputLine )
         
         READ(5,'(A132)') inputLine
         controlVariables % mach = GetRealValue( inputLine )
         
         READ(5,'(A132)') inputLine
         controlVariables % RE = GetRealValue( inputLine )
         
         READ(5,'(A132)') inputLine
         controlVariables % AOATheta = GetRealValue( inputLine )
         
         READ(5,'(A132)') inputLine
         controlVariables % AOAPhi   = GetRealValue( inputLine )
!
!        ---------------------------------------------------------------------------
!        We will store the type and values of the boundaries in dictionaries so that
!        we can associate a name of a boundary curve found in the mesh file with a
!        particular value and type of boundary conditions.
!        ---------------------------------------------------------------------------
!
         READ(5,'(A132)') inputLine
         numberOfBCs = GetIntValue( inputLine )
         
         DO k = 1, numberOfBCs 
            READ(5,*) boundaryName, boundaryValue, boundaryType
            CALL toLower(boundaryName)
            CALL toLower(boundaryType)
            CALL bcTypeDictionary % addValueForKey(boundaryType, boundaryName)
            CALL bcValueDictionary % addValueForKey(boundaryValue, boundaryName)
         END DO
         
      END SUBROUTINE ReadInputFile
