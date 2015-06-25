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
      USE TimeIntegratorClass
      
      IMPLICIT NONE
!
!     ------------
!     Declarations
!     ------------
!
      TYPE( NSLiteControlVariables )      :: controlVariables
      TYPE( DGSem )                       :: sem
      TYPE( FTTimer )                     :: stopWatch
      TYPE( DGSEMPlotter )      , POINTER :: plotter      => NULL()
      CLASS( PlotterDataSource ), POINTER :: plDataSource => NULL()
      TYPE( RKTimeIntegrator )            :: timeIntegrator
      
      REAL(KIND=RP)                       :: dt
      
      LOGICAL                             :: success
      INTEGER                             :: plotUnit, restartUnit
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
      CALL ConstructDGSem(self              = sem, &
                          polynomialOrder   = controlVariables % polynomialOrder,&
                          meshFileName      = controlVariables % inputFileName,  &
                          externalState     = externalStateForBoundaryName,      &
                          externalGradients = ExternalGradientForBoundaryName,   &
                          success           = success)
      IF(.NOT. success)   ERROR STOP "Mesh reading error"
      CALL checkIntegrity(sem % mesh, success)
      IF(.NOT. success)   ERROR STOP "Boundary condition specification error"
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
!     -----------------------------
!     Construct the time integrator
!     -----------------------------
!
      dt = MaxTimeStep( sem, controlVariables % cfl )

      CALL timeIntegrator % constructAsSteadyStateIntegrator(dt            = dt, &
                                                             cfl           = controlVariables % cfl, &
                                                             numberOfSteps = controlVariables % numberOfSteps, &
                                                             plotInterval  = controlVariables % plotInterval)
      CALL timeIntegrator % setIterationTolerance(tol = controlVariables % tol)
!
!     --------------------
!     Prepare for plotting
!     --------------------
!
      IF ( controlVariables % plotFileName /= "none" )     THEN
         plotUnit = UnusedUnit()
         ALLOCATE(plotter)
         ALLOCATE(plDataSource)
         
         CALL plotter % Construct(fUnit      = plotUnit,          &
                                  spA        = sem % spA,         &
                                  dataSource = plDataSource,      &
                                  newN       = controlVariables % numberOfPlotPoints)
         CALL timeIntegrator % setPlotter(plotter)
      END IF 
!
!     -----------------
!     Integrate in time
!     -----------------
!
      CALL stopWatch % start()
      CALL timeIntegrator % integrate(sem)
      CALL stopWatch % stop()
      
      PRINT *, "Elapsed Time: ", stopWatch % elapsedTime(units = TC_SECONDS)
      PRINT *, "Total Time:   ", stopWatch % totalTime(units = TC_SECONDS)
!
!     ----------------
!     Plot the results
!     ----------------
!
      IF ( ASSOCIATED(plotter) )     THEN
         plotUnit = UnusedUnit()
         OPEN(UNIT = plotUnit, FILE = controlVariables % plotFileName)
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
            sem % mesh % elements(eID) % Q(3,3,3,1) = 1.05_RP*sem % mesh % elements(eID) % Q(3,3,3,1)!DEBUG
         END DO 
         
      END SUBROUTINE SetInitialCondition
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE checkIntegrity(mesh, success)
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
               CALL bcObjects % release()
               IF(bcObjects % isUnreferenced()) DEALLOCATE (bcObjects)
               return
            END IF  
         END DO
         
         CALL bcObjects % release()
         IF(bcObjects % isUnreferenced()) DEALLOCATE (bcObjects)
         
      END SUBROUTINE checkIntegrity
