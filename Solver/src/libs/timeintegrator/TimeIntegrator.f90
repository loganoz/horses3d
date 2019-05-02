!
!////////////////////////////////////////////////////////////////////////
!
!      TimeIntegration.f95
!      Created: 2007-10-23 09:25:32 -0400 
!      By: David Kopriva  
!
!      Third order RK integrator for DG approximation to conservation
!      laws in 2D
!
!////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
      MODULE TimeIntegratorClass
      
      USE SMConstants
      use FTValueDictionaryClass
      USE PolynomialInterpAndDerivsModule
      USE DGSEMClass
      use HexMeshClass
      use PhysicsStorage
      USE Physics
      USE ExplicitMethods
      USE IMEXMethods    
      use AutosaveClass                   , only: Autosave_t, AUTOSAVE_BY_TIME
      use StopwatchClass
      use MPI_Process_Info
      use TimeIntegratorDefinitions
      use MonitorsClass
      use ParticlesClass
      use Utilities                       , only: ToLower, AlmostEqual
      use FileReadingUtilities            , only: getFileName
      use ProblemFileFunctions            , only: UserDefinedPeriodicOperation_f
      use pAdaptationClass                , only: pAdaptation_t, ADAPT_DYNAMIC_TIME, ADAPT_STATIC
      use TruncationErrorClass            , only: EstimateAndPlotTruncationError
      use MultiTauEstimationClass         , only: MultiTauEstim_t
      IMPLICIT NONE 
      
      INTEGER, PARAMETER :: TIME_ACCURATE = 0, STEADY_STATE = 1

      TYPE TimeIntegrator_t
         INTEGER                                :: integratorType
         REAL(KIND=RP)                          :: tFinal, time, initial_time
         INTEGER                                :: initial_iter, numTimeSteps, outputInterval, iter
         REAL(KIND=RP)                          :: dt, tolerance, cfl, dcfl
         LOGICAL                                :: Compute_dt                    ! Is st computed from an inputted CFL number?
         type(Autosave_t)                       :: autosave
         type(pAdaptation_t)                    :: pAdaptator
         type(MultiTauEstim_t)                  :: TauEstimator
         PROCEDURE(TimeStep_FCN), NOPASS , POINTER :: RKStep
!
!        ========         
         CONTAINS
!        ========         
!
         PROCEDURE :: construct  => constructTimeIntegrator
         PROCEDURE :: destruct   => destructTimeIntegrator
         PROCEDURE :: integrate
         procedure :: Display    => TimeIntegrator_Display
         procedure :: CorrectDt  => TimeIntegrator_CorrectDt
      END TYPE TimeIntegrator_t

      abstract interface
         subroutine RKStepFcn( sem , t , deltaT )
            use SMConstants, only: RP
            use DGSEMClass,  only: DGSEM
            implicit none
            type(DGSem)     :: sem
            real(kind=RP)   :: t, deltaT
         end subroutine RKStepFcn
      end interface

      character(len=*), parameter   :: TIME_INTEGRATION_KEY  = 'time integration'
      character(len=*), parameter   :: EXPLICIT_SOLVER   = 'explicit'
      character(len=*), parameter   :: IMEX_SOLVER       = 'imex'
      character(len=*), parameter   :: IMPLICIT_SOLVER   = 'implicit'
      character(len=*), parameter   :: FAS_SOLVER        = 'fas'
      character(len=*), parameter   :: ANISFAS_SOLVER    = 'anisfas'
      character(len=*), parameter   :: ROSENBROCK_SOLVER = 'rosenbrock'
!
!     ========      
      CONTAINS 
!     ========      
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE constructTimeIntegrator(self,controlVariables, sem, initial_iter, initial_time)
      
         IMPLICIT NONE
         CLASS(TimeIntegrator_t)     :: self
         TYPE(FTValueDictionary)     :: controlVariables
         type(DGSem)                 :: sem
         integer                     :: initial_iter
         real(kind=RP)               :: initial_time
!
!        ----------------------------------------------------------------------------------
!        Set time-stepping variables
!           If keyword "cfl" is present, the time step size is computed in every time step.
!           If it is not, the keyword "dt" must be specified explicitly.
!        ----------------------------------------------------------------------------------
!
         IF (controlVariables % containsKey("cfl")) THEN
#if defined(NAVIERSTOKES) || defined(INCNS)
            self % Compute_dt = .TRUE.
            self % cfl        = controlVariables % doublePrecisionValueForKey("cfl")
#if defined(NAVIERSTOKES)
            if (flowIsNavierStokes) then
               if (controlVariables % containsKey("dcfl")) then
                  self % dcfl       = controlVariables % doublePrecisionValueForKey("dcfl")
               else
                  ERROR STOP '"cfl" and "dcfl", or "dt" keyword must be specified for the time integrator'
               end if
            end if
#endif
#elif defined(CAHNHILLIARD)
            print*, "Error, use fixed time step to solve Cahn-Hilliard equations"
            errorMessage(STD_OUT)
            stop
#endif
         ELSEIF (controlVariables % containsKey("dt")) THEN
            self % Compute_dt = .FALSE.
            self % dt         = controlVariables % doublePrecisionValueForKey("dt")
         ELSE
            ERROR STOP '"cfl" or "dt" keyword must be specified for the time integrator'
         END IF
!
!        ----------------------
!        Common initializations
!        ----------------------
!
         self % time           =  initial_time 
         self % initial_time   =  initial_time
         self % initial_iter   =  initial_iter
         self % numTimeSteps   =  controlVariables % integerValueForKey ("number of time steps")
         self % outputInterval =  controlVariables % integerValueForKey("output interval")
         self % tolerance      =  controlVariables % doublePrecisionValueForKey("convergence tolerance")
         self % RKStep         => TakeRK3Step
!
!        ------------------------------------
!        Integrator-dependent initializarions
!        ------------------------------------
!
         SELECT CASE (controlVariables % StringValueForKey("simulation type",LINE_LENGTH))
            CASE ('time-accurate')
               IF (controlVariables % containsKey("final time")) THEN
                  self % tFinal         = controlVariables % doublePrecisionValueForKey("final time")
               ELSE
                  self % tFinal         = huge(self % tFinal)
!~                  ERROR STOP '"final time" keyword must be specified for time-accurate integrators'
               ENDIF
               self % integratorType = TIME_ACCURATE
            CASE DEFAULT ! Using 'steady-state' even if not specified in input file
               self % integratorType = STEADY_STATE
         END SELECT
!
!        ----------------
!        Set the autosave
!        ----------------
!
         call self % autosave   % Configure (controlVariables, initial_time)
         call self % pAdaptator % construct (controlVariables, initial_time)      ! If not requested, the constructor returns doing nothing
         
         
         call self % TauEstimator % construct(controlVariables, sem)
         
      END SUBROUTINE constructTimeIntegrator
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE destructTimeIntegrator( self ) 
         CLASS(TimeIntegrator_t) :: self
         self % tFinal       = 0.0_RP
         self % numTimeSteps = 0
         self % dt           = 0.0_RP
         
         if (self % pAdaptator % Constructed) call self % pAdaptator % destruct()
         
         call self % TauEstimator % destruct
      END SUBROUTINE destructTimeIntegrator
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE Integrate( self, sem, controlVariables, monitors, ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
      USE FASMultigridClass
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      CLASS(TimeIntegrator_t)              :: self
      TYPE(DGSem)                          :: sem
      TYPE(FTValueDictionary)              :: controlVariables
      class(Monitor_t)                     :: monitors
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivativeIsolated

!
!     ---------
!     Internal variables
!     ---------
!
      real(kind=RP)        :: FMGres    ! Target residual for FMG solver
      REAL(KIND=RP)        :: maxResidual(NTOTALVARS)
      type(FASMultigrid_t) :: FMGSolver ! FAS multigrid solver for Full-Multigrid (FMG) initialization
      
!     Initializations
!     ---------------

      sem  % numberOfTimeSteps = self % initial_iter
      if (.not. self % Compute_dt) monitors % dt_restriction = DT_FIXED
      
!     Measure solver time
!     -------------------
      
      call Stopwatch % CreateNewEvent("Solver")
      call Stopwatch % Start("Solver")
      
!     Estimate Tau initially, if requested
!     ------------------------------------
      if ( controlVariables % logicalValueForKey("plot truncation error") ) then
         call EstimateAndPlotTruncationError(sem,0._RP,controlVariables,ComputeTimeDerivative,ComputeTimeDerivativeIsolated)
      end if
      
!     Perform FMG cycle if requested
!        (only for steady simulations)
!     ------------------------------
      
      if (controlVariables % containsKey("fasfmg residual")) then
          
         FMGres = controlVariables % doubleprecisionValueForKey("fasfmg residual")
         write(STD_OUT,*) 'Using FMG solver to get initial condition. Res =', FMGres
         
         call FMGSolver % construct(controlVariables,sem)
         call FMGSolver % solve(0,0._RP, 0._RP, ComputeTimeDerivative, .TRUE.,FMGres) 
         
         call FMGSolver % destruct
      end if
      
!     Perform static p-adaptation stage(s) if requested
!     -------------------------------------------------
      if (self % pAdaptator % adaptation_mode == ADAPT_STATIC) then
         
         do while (self % pAdaptator % Adapt)
            
            if (self % integratorType == STEADY_STATE) then
!
!              Lower the residual to 0.1 * truncation error threshold 
!              -> See Kompenhans et al. "Adaptation strategies for high order discontinuous Galerkin methods based on Tau-estimation." Journal of Computational Physics 306 (2016): 216-236.
!              ------------------------------------------------------
               call IntegrateInTime( self, sem, controlVariables, monitors, ComputeTimeDerivative, ComputeTimeDerivativeIsolated, self % pAdaptator % reqTE*0.1_RP)
            end if
            
            call self % pAdaptator % pAdaptTE(sem,sem  % numberOfTimeSteps, self % time, ComputeTimeDerivative, ComputeTimeDerivativeIsolated, controlVariables)
            sem % numberOfTimeSteps = sem % numberOfTimeSteps + 1
            
         end do
      end if
      
!     Finish time integration
!     -----------------------
      call IntegrateInTime( self, sem, controlVariables, monitors, ComputeTimeDerivative, ComputeTimeDerivativeIsolated )
      
!     Measure solver time
!     -------------------
      call Stopwatch % Pause("Solver")

      END SUBROUTINE Integrate    
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------------------
!  Perform the standard time marching integration
!  -> If "tolerance" is provided, the value in controlVariables is ignored. 
!     This is only relevant for STEADY_STATE computations.
!  ------------------------------------------------------------------------
   subroutine IntegrateInTime( self, sem, controlVariables, monitors, ComputeTimeDerivative, ComputeTimeDerivativeIsolated, tolerance, CTD_linear, CTD_nonlinear)
   
      USE BDFTimeIntegrator
      use FASMultigridClass
      use AnisFASMultigridClass
      use RosenbrockTimeIntegrator
      use StopwatchClass
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      CLASS(TimeIntegrator_t)                      :: self
      TYPE(DGSem)                                  :: sem
      TYPE(FTValueDictionary), intent(in)          :: controlVariables
      class(Monitor_t)                             :: monitors
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivativeIsolated
      real(kind=RP), optional, intent(in)          :: tolerance   !< ? tolerance to integrate down to
      procedure(ComputeTimeDerivative_f), optional :: CTD_linear
      procedure(ComputeTimeDerivative_f), optional :: CTD_nonlinear
!
!     ---------------
!     Local variables
!     ---------------
!
      real(kind=RP)                 :: Tol                                 ! Tolerance used for STEADY_STATE computations
      REAL(KIND=RP)                 :: t
      REAL(KIND=RP)                 :: maxResidual(NTOTALVARS)
      REAL(KIND=RP)                 :: dt
      integer                       :: k
      CHARACTER(len=LINE_LENGTH)    :: SolutionFileName
      ! Time-step solvers:
      type(FASMultigrid_t)          :: FASSolver
      type(AnisFASMultigrid_t)      :: AnisFASSolver
      type(BDFIntegrator_t)         :: BDFSolver
      type(RosenbrockIntegrator_t)  :: RosenbrockSolver
      
      CHARACTER(len=LINE_LENGTH)    :: TimeIntegration
      logical                       :: saveGradients
      procedure(UserDefinedPeriodicOperation_f) :: UserDefinedPeriodicOperation
!
!     ----------------------
!     Read Control variables
!     ----------------------
!
      IF (controlVariables % containsKey(TIME_INTEGRATION_KEY)) THEN
         TimeIntegration  = controlVariables % StringValueForKey(TIME_INTEGRATION_KEY,LINE_LENGTH)
      ELSE ! Default value
         TimeIntegration = EXPLICIT_SOLVER
      END IF
      call toLower(TimeIntegration)
      SolutionFileName   = trim(getFileName(controlVariables % StringValueForKey("solution file name",LINE_LENGTH)))
      
!
!     ---------------
!     Initializations
!     ---------------
!
      if (present(tolerance)) then
         Tol = tolerance
      else
         Tol = self % tolerance
      end if
      
      t = self % time
!
!     ------------------
!     Configure restarts
!     ------------------
!
      saveGradients = controlVariables % logicalValueForKey("save gradients with solution")
!
!     -----------------------
!     Check initial residuals
!     -----------------------
!
      call ComputeTimeDerivative(sem % mesh, sem % particles, t, CTD_IGNORE_MODE)
      maxResidual       = ComputeMaxResiduals(sem % mesh)
      sem % maxResidual = maxval(maxResidual)
      call Monitors % UpdateValues( sem % mesh, t, sem % numberOfTimeSteps, maxResidual )
      call self % Display(sem % mesh, monitors, sem  % numberOfTimeSteps)
      
      if (self % pAdaptator % adaptation_mode    == ADAPT_DYNAMIC_TIME .and. &
          self % pAdaptator % nextAdaptationTime == self % time) then
         call self % pAdaptator % pAdaptTE(sem,sem  % numberOfTimeSteps,t, ComputeTimeDerivative, ComputeTimeDerivativeIsolated, controlVariables)
         self % pAdaptator % nextAdaptationTime = self % pAdaptator % nextAdaptationTime + self % pAdaptator % time_interval
      end if 
      
      call monitors % WriteToFile(sem % mesh)

      IF (self % integratorType == STEADY_STATE) THEN
         IF (maxval(maxResidual) <= Tol )  THEN
            if (MPI_Process % isRoot) then
               write(STD_OUT,'(/,A,I0,A,ES10.3)') "   *** Residual tolerance reached at iteration ",sem % numberOfTimeSteps," with Residual = ", maxval(maxResidual)
            end if
            call monitors % WriteToFile(sem % mesh, force = .TRUE.)
            return
         END IF
      end if
!
!     -----------------
!     Integrate in time
!     -----------------
!
      select case (TimeIntegration)
      case(FAS_SOLVER)
         call FASSolver % construct(controlVariables,sem)

      case(ANISFAS_SOLVER)
         call AnisFASSolver % construct(controlVariables,sem)

      case(IMPLICIT_SOLVER)
         call BDFSolver % construct(controlVariables,sem)

      case(ROSENBROCK_SOLVER)
         call RosenbrockSolver % construct(controlVariables,sem)

      end select
      
      DO k = sem  % numberOfTimeSteps, self % initial_iter + self % numTimeSteps-1
!
!        CFL-bounded time step
!        ---------------------      
         IF ( self % Compute_dt ) self % dt = MaxTimeStep( sem, self % cfl, self % dcfl )
!
!        Correct time step
!        -----------------
         dt = self % CorrectDt(t,self % dt)
!
!        User defined periodic operation
!        -------------------------------
         CALL UserDefinedPeriodicOperation(sem % mesh, t, dt, monitors)
!
!        Perform time step
!        -----------------         
         SELECT CASE (TimeIntegration)
         CASE (IMPLICIT_SOLVER)
            call BDFSolver % TakeStep (sem, t , dt , ComputeTimeDerivative)
         CASE (ROSENBROCK_SOLVER)
            call RosenbrockSolver % TakeStep (sem, t , dt , ComputeTimeDerivative)
         CASE (EXPLICIT_SOLVER)
            CALL self % RKStep ( sem % mesh, sem % particles, t, dt, ComputeTimeDerivative)
         case (FAS_SOLVER)
            call FASSolver % solve(k, t, dt, ComputeTimeDerivative)
         case (ANISFAS_SOLVER)
            call AnisFASSolver % solve(k,t, ComputeTimeDerivative)
         case (IMEX_SOLVER)
            call TakeIMEXStep(sem, t, dt, controlVariables, computeTimeDerivative)
         END SELECT
!
!        Compute the new time
!        --------------------         
         t = t + dt
         self % time = t
!
!        Get maximum residuals
!        ---------------------
         maxResidual       = ComputeMaxResiduals(sem % mesh)
         sem % maxResidual = maxval(maxResidual)
!
!        Update monitors
!        ---------------
         call Monitors % UpdateValues( sem % mesh, t, k+1, maxResidual )
!
!        Exit if the target is reached
!        -----------------------------
         IF (self % integratorType == STEADY_STATE) THEN
            IF (maxval(maxResidual) <= Tol )  THEN
               call self % Display(sem % mesh, monitors, k+1)
               if (MPI_Process % isRoot) then
                  write(STD_OUT,'(/,A,I0,A,ES10.3)') "   *** Residual tolerance reached at iteration ",k+1," with Residual = ", maxval(maxResidual)
               end if
               sem % numberOfTimeSteps = k + 1
               exit
            END IF
         ELSEIF (self % integratorType == TIME_ACCURATE) THEN
            IF ( t .ge. self % tFinal) then
               call self % Display( sem % mesh, monitors, k+1)
               sem % numberOfTimeSteps = k + 1
               exit
            end if
         END IF
#if defined(NAVIERSTOKES)
!
!        Integration of particles
!        ------------------------
         if ( sem % particles % active ) then 
            call sem % particles % Integrate(sem % mesh, dt)
         endif 
#endif
!
!        Print monitors
!        --------------
         IF( (MOD( k+1, self % outputInterval) == 0) .or. (k .eq. self % initial_iter) ) call self % Display(sem % mesh, monitors, k+1)
!
!        p- Adapt
!        --------------
         IF( self % pAdaptator % hasToAdapt(k+1) ) then
            call self % pAdaptator % pAdaptTE(sem,k,t, ComputeTimeDerivative, ComputeTimeDerivativeIsolated, controlVariables)
         end if
         call self % TauEstimator % estimate(sem, k+1, t, ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
!
!        Autosave
!        --------         
         if ( self % autosave % Autosave(k+1) ) then
            call SaveRestart(sem,k+1,t,SolutionFileName, saveGradients)
   
         end if

!        Flush monitors
!        --------------
         call monitors % WriteToFile(sem % mesh)
         
         sem % numberOfTimeSteps = k + 1
      END DO
!
!     Flush the remaining information in the monitors
!     -----------------------------------------------
      if ( k .ne. 0 ) then
         call Monitors % writeToFile(sem % mesh, force = .true. )
      end if
      
      sem % maxResidual       = maxval(maxResidual)
      self % time             = t
      
!
!     ---------
!     Finish up
!     ---------
!
      select case(TimeIntegration)
      case(FAS_SOLVER)
         CALL FASSolver % destruct
      
      case(ANISFAS_SOLVER)
         CALL AnisFASSolver % destruct
      
      case(IMPLICIT_SOLVER)
         call BDFSolver % destruct

      case(ROSENBROCK_SOLVER)
         call RosenbrockSolver % destruct

      end select

   end subroutine IntegrateInTime
      
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
!     Subroutine to print the residuals
!
!
   subroutine TimeIntegrator_Display(self, mesh, monitors, iter)
      implicit none
      class(TimeIntegrator_t),   intent(in)     :: self
      class(HexMesh),            intent(in)     :: mesh
      class(Monitor_t),          intent(inout)  :: monitors
      integer                  , intent(in)     :: iter
!
!     ---------------
!     Local variables      
!     ---------------
!
      real(kind=RP)           :: ETA, tEl
      integer, parameter      :: showLabels = 50
      integer, save           :: shown = 0

      if ( .not. MPI_Process % isRoot ) return 

      if ( mod(shown, showLabels) .eq. 0 ) then
         if ( (self % integratorType .eq. TIME_ACCURATE) .and. (iter .gt. self % initial_iter+1) ) then 
!
!           Compute ETA
!           -----------
            tEl = Stopwatch % ElapsedTime("Solver")
            ETA = (self % tFinal - self % initial_time) * tEl / (self % time - self % initial_time) - tEl
            write(STD_OUT,'(A,F10.3,A)') "*** ETA:", ETA," seconds."
         end if
         write(STD_OUT,'(/)')
         write(STD_OUT,'(/)')
         
         call monitors % WriteLabel
         call monitors % WriteUnderlines

      end if
      shown = shown + 1 

      call monitors % WriteValues

   end subroutine TimeIntegrator_Display
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE SaveRestart(sem,k,t,RestFileName, saveGradients)
      IMPLICIT NONE
!
!     ------------------------------------
!     Save the results to the restart file
!     ------------------------------------
!
!     ----------------------------------------------
      TYPE(DGSem)                  :: sem            !< DGsem class
      INTEGER                      :: k              !< Time step
      REAL(KIND=RP)                :: t              !< Simu time
      CHARACTER(len=*)             :: RestFileName   !< Name of restart file
      logical,          intent(in) :: saveGradients
!     ----------------------------------------------
      INTEGER                      :: fd             !  File unit for new restart file
      CHARACTER(len=LINE_LENGTH)   :: FinalName      !  Final name for particular restart file
!     ----------------------------------------------
      
      WRITE(FinalName,'(2A,I10.10,A)')  TRIM(RestFileName),'_',k,'.hsol'
      if ( MPI_Process % isRoot ) write(STD_OUT,'(A,A,A,ES10.3,A)') '*** Writing file "',trim(FinalName),'", with t = ',t,'.'
      call sem % mesh % SaveSolution(k,t,trim(finalName),saveGradients)
   
   END SUBROUTINE SaveRestart
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  This routine corrects the time-step size, so that 
!  time-periodic operations can be performed
!  -------------------------------------------------
   recursive function TimeIntegrator_CorrectDt (self, t, dt_in) result(dt_out)
      implicit none
      !-arguments------------------------------------------------
      class(TimeIntegrator_t) , intent(inout)   :: self
      real(kind=RP)           , intent(in)      :: t
      real(kind=RP)           , intent(in)      :: dt_in
      real(kind=RP)                             :: dt_out
      !-local-variables------------------------------------------
      real(kind=RP) :: dt_temp
      
      integer, parameter :: DO_NOTHING  = 0
      integer, parameter :: AUTOSAVE    = 1
      integer, parameter :: ADAPT       = 2
      integer, parameter :: DONT_KNOW   = 3
      integer, save :: next_time_will = DONT_KNOW
      !----------------------------------------------------------

!
!     Initializations
!     -------------------------------      
      self % pAdaptator % performPAdaptationT = .FALSE.
      self % autosave   % performAutosave = .FALSE.
      dt_out = dt_in
      
!
!     time-step bounded by final time
!     -------------------------------
      if ( self % integratorType .eq. TIME_ACCURATE ) then
         if ( ( t + dt_out) .gt. self % tFinal ) then
            dt_out = self % tFinal - t
         end if
      end if
      
!
!     time-step bounded by periodic operations
!     ----------------------------------------
      
      select case (next_time_will)
         case (DO_NOTHING)
            return
            
         case (AUTOSAVE)
            
            if ( self % autosave % nextAutosaveTime < (t + dt_out) ) then
               dt_out = self % autosave % nextAutosaveTime - t
               self % autosave % performAutosave = .TRUE.
               
               if ( AlmostEqual(self % autosave % nextAutosaveTime, self % pAdaptator % nextAdaptationTime) ) then
                  self % pAdaptator % performPAdaptationT = .TRUE.
                  self % pAdaptator % nextAdaptationTime = self % pAdaptator % nextAdaptationTime + self % pAdaptator % time_interval
               end if
               
               self % autosave % nextAutosaveTime = self % autosave % nextAutosaveTime + self % autosave % time_interval
               next_time_will = minloc([self % autosave % nextAutosaveTime, self % pAdaptator % nextAdaptationTime],1)
            end if
            
         case (ADAPT)
            
            if ( self % pAdaptator % nextAdaptationTime < (t + dt_out) ) then
               dt_out = self % pAdaptator % nextAdaptationTime - t
               self % pAdaptator % performPAdaptationT = .TRUE.
               
               if ( AlmostEqual(self % autosave % nextAutosaveTime, self % pAdaptator % nextAdaptationTime) ) then
                  self % autosave % performAutosave = .TRUE.
                  self % autosave % nextAutosaveTime = self % autosave % nextAutosaveTime + self % autosave % time_interval
               end if
               
               self % pAdaptator % nextAdaptationTime = self % pAdaptator % nextAdaptationTime + self % pAdaptator % time_interval
               next_time_will = minloc([self % autosave % nextAutosaveTime, self % pAdaptator % nextAdaptationTime],1)
            end if
            
         case (DONT_KNOW)
            
            if (  self % pAdaptator % adaptation_mode == ADAPT_DYNAMIC_TIME .or. &
                  self % autosave % mode       == AUTOSAVE_BY_TIME) then
               
               next_time_will = minloc([self % autosave % nextAutosaveTime, self % pAdaptator % nextAdaptationTime],1)
               
               dt_temp = self % CorrectDt (t, dt_out)
               dt_out  = dt_temp
            else
               next_time_will = DO_NOTHING
            end if
            
      end select
      
   end function TimeIntegrator_CorrectDt
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!         
END MODULE TimeIntegratorClass
