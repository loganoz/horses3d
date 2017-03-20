!
!////////////////////////////////////////////////////////////////////////
!
!      NSLite3D.f90
!      Created: May 21, 2015 at 12:56 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
!
!////////////////////////////////////////////////////////////////////////
!
      PROGRAM NSLite3DMain
      
      USE SMConstants
      USE FTTimerClass
      USE PhysicsStorage
      USE SharedBCModule
      USE DGSEMPlotterClass
      USE DGSEMClass
      USE BoundaryConditionFunctions
      USE TimeIntegratorClass
      USE UserDefinedFunctions
      USE mainKeywordsModule
      USE Headers
      
      IMPLICIT NONE
!
!     ------------
!     Declarations
!     ------------
!
      TYPE( FTValueDictionary)            :: controlVariables
      TYPE( DGSem )                       :: sem
      TYPE( FTTimer )                     :: stopWatch
      TYPE( DGSEMPlotter )      , POINTER :: plotter      => NULL()
      CLASS( PlotterDataSource ), POINTER :: plDataSource => NULL()
      TYPE( TimeIntegrator_t )            :: timeIntegrator
      
      REAL(KIND=RP)                       :: dt, cfl
      LOGICAL                             :: success
      INTEGER                             :: plotUnit, restartUnit, saveUnit
      INTEGER, EXTERNAL                   :: UnusedUnit
      EXTERNAL                            :: externalStateForBoundaryName
      EXTERNAL                            :: ExternalGradientForBoundaryName
!
!     ---------------
!     Initializations
!     ---------------
!
      CALL Main_Header("NSLite3D discontinuous Galerkin spectral element solver.")

      CALL controlVariables % initWithSize(16)
      CALL stopWatch % init()
      CALL UserDefinedStartup
      CALL ConstructSharedBCModule
      
      CALL ReadInputFile( controlVariables )
      CALL CheckInputIntegrity(controlVariables, success)
      IF(.NOT. success)   ERROR STOP "Control file reading error"
      
!
!     ----------------
!     Set up the DGSEM
!     ----------------
!      
      CALL ConstructPhysicsStorage( controlVariables, success )
      IF(.NOT. success)   ERROR STOP "Physics parameters input error"
                                   
      CALL sem % construct(polynomialOrder   = controlVariables % integerValueForKey(polynomialOrderKey),&
                           meshFileName      = controlVariables % stringValueForKey(meshFileNameKey,     &
                                                                        requestedLength = LINE_LENGTH),  &
                           externalState     = externalStateForBoundaryName,                             &
                           externalGradients = ExternalGradientForBoundaryName,                          &
                           success           = success)
                           
      IF(.NOT. success)   ERROR STOP "Mesh reading error"
      CALL checkBCIntegrity(sem % mesh, success)
      IF(.NOT. success)   ERROR STOP "Boundary condition specification error"
      CALL UserDefinedFinalSetup(sem, controlVariables)
!
!     ----------------------
!     Set the initial values
!     ----------------------
!
      IF ( controlVariables % logicalValueForKey(restartKey) )     THEN
         restartUnit = UnusedUnit()
         OPEN( UNIT = restartUnit, &
               FILE = controlVariables % stringValueForKey(restartFileNameKey,requestedLength = LINE_LENGTH), &
               FORM = "UNFORMATTED" )
               CALL sem % LoadSolutionForRestart( restartUnit )
         CLOSE( restartUnit )
      ELSE
         CALL UserDefinedInitialCondition(sem , controlVariables)
      END IF
!
!     -----------------------------
!     Construct the time integrator
!     -----------------------------
!
     CALL timeIntegrator % construct (sem,controlVariables)
!
!     --------------------
!     Prepare for plotting
!     --------------------
!
      IF ( controlVariables % stringValueForKey(plotFileNameKey, &
           requestedLength = LINE_LENGTH) /= "none" )     THEN
         plotUnit = UnusedUnit()
         ALLOCATE(plotter)
         ALLOCATE(plDataSource)
         
         CALL plotter % Construct(fUnit      = plotUnit,          &
                                  spA        = sem % spA,         &
                                  dataSource = plDataSource,      &
                                  newN       = controlVariables % integerValueForKey(numberOfPlotPointsKey))
         CALL timeIntegrator % setPlotter(plotter)
      END IF 
!
!     -----------------
!     Integrate in time
!     -----------------
!
      CALL stopWatch % start()
         CALL timeIntegrator % integrate(sem, controlVariables)
      CALL stopWatch % stop()
      
      PRINT *
      PRINT *, "Elapsed Time: ", stopWatch % elapsedTime(units = TC_SECONDS)
      PRINT *, "Total Time:   ", stopWatch % totalTime  (units = TC_SECONDS)
!
!     -----------------------------------------------------
!     Let the user perform actions on the computed solution
!     -----------------------------------------------------
!
      CALL UserDefinedFinalize(sem, timeIntegrator % time)
!
!     ------------------------------------
!     Save the results to the restart file
!     ------------------------------------
!
      IF(controlVariables % stringValueForKey(saveFileNameKey,LINE_LENGTH) /= "none")     THEN 
         saveUnit = UnusedUnit()
         OPEN( UNIT = saveUnit, &
               FILE = controlVariables % stringValueForKey(saveFileNameKey,LINE_LENGTH), &
               FORM = "UNFORMATTED" )
               CALL sem % SaveSolutionForRestart( saveUnit )
         CLOSE( saveUnit )
      END IF
!
!     ----------------
!     Plot the results
!     ----------------
!
      IF ( ASSOCIATED(plotter) )     THEN
         plotUnit = UnusedUnit()
         OPEN(UNIT = plotUnit, FILE = controlVariables % stringValueForKey(plotFileNameKey, &
                                                                requestedLength = LINE_LENGTH))
            CALL plotter % ExportToTecplot( elements = sem % mesh % elements )
         CLOSE(plotUnit)
      END IF 
!
!     --------
!     Clean up
!     --------
!
      IF(ASSOCIATED(plotter)) THEN
         CALL plotter % Destruct()
         DEALLOCATE(plotter)
         DEALLOCATE(plDataSource)
      END IF 
      CALL timeIntegrator % destruct()
      CALL sem % destruct()
      CALL destructSharedBCModule
      
      CALL UserDefinedTermination
      
      END PROGRAM NSLite3DMain
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE CheckBCIntegrity(mesh, success)
!
         USE HexMeshClass
         USE SharedBCModule
         USE BoundaryConditionFunctions, ONLY:implementedBCNames
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh) :: mesh
         LOGICAL       :: success
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                              :: i, j
         INTEGER                              :: faceID, eId
         CHARACTER(LEN=BC_STRING_LENGTH)      :: bcName, namedBC
         CHARACTER(LEN=BC_STRING_LENGTH)      :: bcType
         CLASS(FTMutableObjectArray), POINTER :: bcObjects
         CLASS(FTValue)             , POINTER :: v
         CLASS(FTObject), POINTER             :: obj
         
         success = .TRUE.
!
!        ----------------------------------------------------------
!        Check to make sure that the boundaries defined in the mesh
!        have an associated name in the control file.
!        ----------------------------------------------------------
         
         DO eID = 1, SIZE( mesh % elements )
            DO faceID = 1, 6
               namedBC = mesh % elements(eId) % boundaryName(faceID)
               IF( namedBC == emptyBCName ) CYCLE
               
               bcName = bcTypeDictionary % stringValueForKey(key             = namedBC,         &
                                                             requestedLength = BC_STRING_LENGTH)
               IF ( LEN_TRIM(bcName) == 0 )     THEN
                  PRINT *, "Control file does not define a boundary condition for boundary name = ", &
                            mesh % elements(eId) % boundaryName(faceID)
                  success = .FALSE.
                  return 
               END IF 
            END DO   
         END DO
!
!        --------------------------------------------------------------------------
!        Check that the boundary conditions to be applied are implemented
!        in the code. Keep those updated in the boundary condition functions module
!        --------------------------------------------------------------------------
!
         bcObjects => bcTypeDictionary % allObjects()
         DO j = 1, bcObjects % COUNT()
            obj => bcObjects % objectAtIndex(j)
            CALL castToValue(obj,v)
            bcType = v % stringValue(requestedLength = BC_STRING_LENGTH)
            DO i = 1, SIZE(implementedBCNames)
               IF ( bcType == implementedBCNames(i) )     THEN
                  success = .TRUE. 
                  EXIT 
               ELSE 
                  success = .FALSE. 
               END IF 
            END DO
            
            IF ( .NOT. success )     THEN
               PRINT *, "Boundary condition ", TRIM(bcType)," not implemented in this code"
               CALL release(bcObjects)
               RETURN 
            END IF  
            
         END DO
         
         CALL release(bcObjects)
         
      END SUBROUTINE checkBCIntegrity
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE CheckInputIntegrity( controlVariables, success )  
         USE FTValueDictionaryClass
         USE mainKeywordsModule
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(FTValueDictionary) :: controlVariables
         LOGICAL                 :: success
!
!        ---------------
!        Local variables
!        ---------------
!
         CLASS(FTObject), POINTER :: obj
         INTEGER                  :: i
         success = .TRUE.
         
         DO i = 1, SIZE(mainKeywords)
            obj => controlVariables % objectForKey(mainKeywords(i))
            IF ( .NOT. ASSOCIATED(obj) )     THEN
               PRINT *, "Input file is missing entry for keyword: ",mainKeywords(i)
               success = .FALSE. 
            END IF  
         END DO  
         
         
      END SUBROUTINE checkInputIntegrity
