!
!////////////////////////////////////////////////////////////////////////
!
!      NSLite3D.f90
!      Created: May 21, 2015 at 12:56 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      PROGRAM NSLite3DMain
      
      USE SMConstants
      USE FTTimerClass
      USE PhysicsStorage
      USE SharedBCModule
      USE ControlVariablesModule
      USE DGSEMPlotterClass
      USE DGSEMClass
      USE BoundaryConditionFunctions
      
      IMPLICIT NONE
!
!     ------------
!     Declarations
!     ------------
!
      TYPE( NSLiteControlVariables )      :: controlVariables
      TYPE( DGSem )                       :: sem
      TYPE( FTTimer )                     :: stopWatch
      TYPE( DGSEMPlotter )                :: plotter
      CLASS( PlotterDataSource ), POINTER :: plDataSource
      
      LOGICAL                             :: success
      INTEGER                             :: plotUnit, restartUnit
      INTEGER                             :: numberOfPlotPoints = 5
      INTEGER, EXTERNAL                   :: UnusedUnit
      EXTERNAL                            :: externalStateForBoundaryName
      EXTERNAL                            :: ExternalGradientForBoundaryName
!
!     ---------------
!     Initializations
!     ---------------
!
      CALL constructSharedBCModule
      CALL ReadInputFile( controlVariables )
      CALL ConstructPhysicsStorage( controlVariables % mach,               &
                                    controlVariables % RE,                 &
                                    0.72_RP,                               &
                                    controlVariables % flowIsNavierStokes )
      CALL stopWatch % init()
!
!     ----------------
!     Set up the DGSEM
!     ----------------
!
      CALL stopWatch % start()
      
      CALL ConstructDGSem(self              = sem, &
                          polynomialOrder   = controlVariables % polynomialOrder,&
                          meshFileName      = controlVariables % inputFileName,  &
                          externalState     = externalStateForBoundaryName,      &
                          externalGradients = ExternalGradientForBoundaryName,   &
                          success           = success)
!
!     ----------------------
!     Set the initial values
!     ----------------------
!
      IF ( controlVariables % restart )     THEN
         restartUnit = UnusedUnit()
         OPEN( UNIT = restartUnit, &
               FILE = controlVariables % restartFileName, &
               FORM = "UNFORMATTED" )
               CALL LoadSolutionForRestart( sem, restartUnit )
         CLOSE( restartUnit )
      ELSE
         CALL SetInitialCondition( sem, UniformFlowState )
      END IF
!
!     ----------------
!     Plot the results
!     ----------------
!
      plotUnit = UnusedUnit()
      ALLOCATE(plDataSource)
      
      OPEN(UNIT=plotUnit, FILE = controlVariables % plotFileName)
         CALL plotter % Construct(fUnit      = plotUnit,          &
                                  spA        = sem % spA,         &
                                  dataSource = plDataSource,      &
                                  newN       = numberOfPlotPoints)
         CALL plotter % ExportToTecplot( elements = sem % mesh % elements )
      CLOSE(plotUnit)
      
      CALL stopWatch % stop()
      PRINT *, stopWatch % elapsedTime(units = TC_SECONDS), &
               stopWatch % totalTime(units = TC_SECONDS)
!
!     --------
!     Clean up
!     --------
!
      CALL plotter % Destruct()
      DEALLOCATE(plDataSource)
      CALL sem % destruct()
      CALL destructSharedBCModule
      
      END PROGRAM NSLite3DMain
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SetInitialCondition( sem, initialStateSubroutine )
         USE SMConstants
         USE DGSEMClass
         USE PhysicsStorage
         IMPLICIT NONE
         
         TYPE(DGSem)      :: sem
         EXTERNAL         :: initialStateSubroutine
                  
         INTEGER     :: i, j, k, eID
         
         DO eID = 1, SIZE(sem % mesh % elements)
            DO k = 0, sem % spA % N
               DO j = 0, sem % spA % N
                  DO i = 0, sem % spA % N 
                     CALL initialStateSubroutine( sem % mesh % elements(eID) % geom % x(:,i,j,k), 0.0_RP, &
                                                  sem % mesh % elements(eID) % Q(i,j,k,1:N_EQN) )
                  END DO
               END DO
            END DO 
         END DO 
         
      END SUBROUTINE SetInitialCondition
