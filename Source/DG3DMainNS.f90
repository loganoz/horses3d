!
!////////////////////////////////////////////////////////////////////////
!
!      DG3DMainNS.f90
!      Created: 2007-10-22 11:36:38 -0400 
!      By: David Kopriva

!>     Compute the solution of a conservation law on a mapped domain using
!!      a discontinuous Galerkin spectral element approxmation. 
!
!////////////////////////////////////////////////////////////////////////
!
      Module ControlVariablesModule
         USE SMConstants
         IMPLICIT NONE
         
         CHARACTER(LEN=LINE_LENGTH) :: inputFileName, plotFileName, restartFileName
         INTEGER                    :: plotInterval
         INTEGER                    :: numberOfSteps
         LOGICAL                    :: restart
         
      END MODULE ControlVariablesModule
!
!////////////////////////////////////////////////////////////////////////
!
!                                   MAIN
!
!////////////////////////////////////////////////////////////////////////
!      
      PROGRAM DGSemMain
!
!     ------------
!     Modules used
!     ------------
!
      USE SMConstants
      USE DGSEMClass
      USE TimeIntegratorClass
      USE PDEModule
      USE TransfiniteMapClass
      USE SharedBCModule
      USE BoundaryConditionFunctions
      USE PlotterClass
      USE ControlVariablesModule
      USE FTTimerClass
      
      IMPLICIT NONE 
!
      TYPE(TimeIntegrator) :: integrator
      TYPE(DGSem)          :: sem
      TYPE(Plotter)        :: thePlotter
      
      INTEGER              :: N, fUnit = 11
      REAL(KIND=RP)        :: dt, cfl
      REAL(KIND=RP)        :: tol
      
      CHARACTER(LEN=1)     :: tab = CHAR(9)
!
!
!     ------
!     Timing
!     ------
!
      TYPE( FTTimer ) :: stopWatch
!
!     ------------------
!     External functions
!     ------------------
!
      EXTERNAL                :: ExternalStateForBoundaryName, ExternalGradientForBoundaryName
      REAL(KIND=RP), EXTERNAL :: MaximumEigenvalue, EstimateMaximumEigenvalue
!
      CALL stopWatch  %  init()
!
      CALL ReadInputFile(N, cfl, tol)
      CALL ConstructPhysicsStorage( mach, RE, 0.72_RP)
      riemannSolverChoice = ROE
!
!     ---------------------------
!     Setup spatial approximation
!     and initial condition.
!     ---------------------------
!     
      CALL ConstructDGSem( sem, N, inputFileName )
      IF ( restart )     THEN
         OPEN( UNIT = 11, FILE = restartFileName, FORM = "UNFORMATTED" )
            CALL LoadSolutionForRestart( sem, 11 )
         CLOSE( 11 )
      ELSE
         CALL SetInitialCondition( sem, UniformFlowState )
      END IF
      
!
!     -------------------------------------------------------------
!     Set up time integrator. With these arguments, the time
!     integrator is being set up to compute a steady-state problem.
!     See the time integrator class source.
!     -------------------------------------------------------------
!
      dt         = MaxTimeStep( sem, cfl )
!
      integrator = NewTimeIntegrator( dt, numberOfSteps, plotInterval )
      CALL SetIterationTolerance( integrator, tol )
!
!     -------------------
!     Set up for plotting
!     -------------------
!
      CALL ConstructInterpolatingPlotter( thePlotter, fUnit, 12, 12 )
!
!     -----------------
!     Integrate in time
!     -----------------
!
      CALL stopWatch  %  start()
         CALL Integrate( integrator, sem, cfl, thePlotter, &
                         ExternalStateForBoundaryName, ExternalGradientForBoundaryName )
      CALL stopWatch  %  stop()
      PRINT *, "polynomial Order = ", N, CHAR(9), " # of steps = ", &
               numberOfSteps, CHAR(9), "Wall Clock Time (min.) = ", stopWatch  %  elapsedTime(TC_MINUTES), &
               "Total CPU Time(min.) = ", stopWatch  %  totalTime(TC_MINUTES)
!
!     ---------
!     Finish up
!     ---------
!
      CALL Destruct( integrator )
!
!     --------------------------------
!     Save results in the restart file
!     --------------------------------
!
      OPEN(UNIT = fUnit, FILE = restartFileName, FORM = "UNFORMATTED" )
         CALL SaveSolutionForRestart( sem, fUnit )
      CLOSE( fUnit )
!
!     --------------
!     Output results
!     --------------
!
      OPEN( UNIT= fUnit, FILE = plotFileName )
         CALL ExportToTecplot( thePlotter, sem )
      CLOSE(fUnit)
!
!     --------------
!     Compute errors
!     --------------
!
!      CALL ComputeErrors(sem)
      CALL Destruct( sem )
      
      END PROGRAM DGSemMain
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SetInitialCondition( sem, solution )
         USE SMConstants
         USE DGSEMClass
         USE PDEModule
         IMPLICIT NONE
         
         TYPE(DGSem)      :: sem
         EXTERNAL         :: solution
                  
         INTEGER     :: i, j, eID, k
         
         k = sem % spA % N/2
         
         DO eID = 1, SIZE(sem % mesh % elements)
            DO j = 0, sem % spA % N
               DO i = 0, sem % spA % N 
                  CALL solution( sem % mesh % elements(eID) % geom % x(i,j), &
                                 sem % mesh % elements(eID) % geom % y(i,j), 0.0_RP, &
                                 sem % mesh % elements(eID) % Q(i,j,1:nEqn) )
               END DO
            END DO
         END DO 
         
      END SUBROUTINE SetInitialCondition
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE ReadInputFile(polynomialOrder,cfl, tol)
         USE SharedBCModule  
         USE PDEModule
         USE ControlVariablesModule
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         INTEGER       :: polynomialOrder
         REAL(KIND=RP) :: cfl
         REAL(KIND=RP) :: tol
!
!        ---------------
!        Local variables
!        ---------------
!
         CHARACTER(LEN=LINE_LENGTH)              :: inputLine
         INTEGER                                 :: numberOfBCs, k
         CHARACTER(LEN=DICT_KWD_STRING_LENGTH)   :: boundaryName
         CHARACTER(LEN=DICT_KWD_STRING_LENGTH)   :: flowEquationsName
         CHARACTER(LEN=DICT_KWD_STRING_LENGTH)   :: equationformName
         CHARACTER(LEN=DICT_VALUE_STRING_LENGTH) :: boundaryType
         CHARACTER(LEN=DICT_VALUE_STRING_LENGTH) :: boundaryValue
!
!        ------------------
!        External functions
!        ------------------
!
         REAL(KIND=RP)             , EXTERNAL    :: GetRealValue
         INTEGER                   , EXTERNAL    :: GetIntValue
         CHARACTER(LEN=LINE_LENGTH), EXTERNAL    :: GetStringValue
         LOGICAL                   , EXTERNAL    :: GetLogicalValue

         CALL ConstructDictionary( bcValueDictionary )
         CALL ConstructDictionary( bcTypeDictionary )
!
!        -----------------------------------------------
!        Read the input file. 
!        Nothing fancy here. Pretty much a fixed format.
!        Except that for convenience, we use a dictionary
!        to keep what BC type is to be applied to which 
!        curve name.
!        -----------------------------------------------
!
         READ(5,'(A132)') inputLine
         flowEquationsName = GetStringValue( inputLine )
         
         IF ( flowEquationsName == 'Euler' .OR. flowEquationsName == 'euler' )     THEN
            flowIsNavierStokes = .false.
         ELSE
            flowIsNavierStokes = .true.
         END IF
   
         READ(5,'(A132)') inputLine
         inputFileName = GetStringValue( inputLine )
         
         READ(5,'(A132)') inputLine
         plotFileName = GetStringValue( inputLine )
         
         READ(5,'(A132)') inputLine
         restartFileName = GetStringValue( inputLine )
   
         READ(5,'(A132)') inputLine
         restart = GetLogicalValue( inputLine )
         
         READ(5,'(A132)') inputLine
         polynomialOrder = GetIntValue( inputLine )
         
         READ(5,'(A132)') inputLine
         numberOfSteps = GetIntValue( inputLine )
         
         READ(5,'(A132)') inputLine
         plotInterval = GetIntValue( inputLine )
         
         READ(5,'(A132)') inputLine
         tol = GetRealValue( inputLine )
         
         READ(5,'(A132)') inputLine
         cfl = GetRealValue( inputLine )
         
         READ(5,'(A132)') inputLine
         mach = GetRealValue( inputLine )
         
         READ(5,'(A132)') inputLine
         RE = GetRealValue( inputLine )
         
         READ(5,'(A132)') inputLine
         AOA = GetRealValue( inputLine )         
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
            CALL AddValue_ForKey_ToDict_( boundaryType , boundaryName, bcTypeDictionary )
            CALL AddValue_ForKey_ToDict_( boundaryValue, boundaryName, bcValueDictionary )
         END DO
      END SUBROUTINE ReadInputFile
