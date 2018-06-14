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
      use AutosaveClass
      use StopwatchClass
      use MPI_Process_Info
      use TimeIntegratorDefinitions
      use MonitorsClass
      use ParticlesClass
      use Utilities, only: ToLower
      use FileReadingUtilities      , only: getFileName
      IMPLICIT NONE 
      
      INTEGER, PARAMETER :: TIME_ACCURATE = 0, STEADY_STATE = 1

      TYPE TimeIntegrator_t
         INTEGER                                :: integratorType
         REAL(KIND=RP)                          :: tFinal, time, initial_time
         INTEGER                                :: initial_iter, numTimeSteps, outputInterval, iter
         REAL(KIND=RP)                          :: dt, tolerance, cfl, dcfl
         LOGICAL                                :: Compute_dt                    ! Is st computed from an inputted CFL number?
         type(Autosave_t)                       :: autosave
         PROCEDURE(TimeStep_FCN), NOPASS , POINTER :: RKStep
!
!        ========         
         CONTAINS
!        ========         
!
         PROCEDURE :: construct => constructTimeIntegrator
         PROCEDURE :: destruct => destructTimeIntegrator
         PROCEDURE :: integrate
         procedure :: Display => TimeIntegrator_Display
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
      SUBROUTINE constructTimeIntegrator(self,controlVariables, initial_iter, initial_time)
      
         IMPLICIT NONE
         CLASS(TimeIntegrator_t)     :: self
         TYPE(FTValueDictionary)     :: controlVariables
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
#if defined(NAVIERSTOKES)
            self % Compute_dt = .TRUE.
            self % cfl        = controlVariables % doublePrecisionValueForKey("cfl")
            if (flowIsNavierStokes) then
               if (controlVariables % containsKey("dcfl")) then
                  self % dcfl       = controlVariables % doublePrecisionValueForKey("dcfl")
               else
                  ERROR STOP '"cfl" and "dcfl", or "dt" keyword must be specified for the time integrator'
               end if
            end if
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
                  ERROR STOP '"final time" keyword must be specified for time-accurate integrators'
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
         call self % autosave % Configure(controlVariables, initial_time)
         
      END SUBROUTINE constructTimeIntegrator
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE destructTimeIntegrator( self ) 
         CLASS(TimeIntegrator_t) :: self
         self % tFinal       = 0.0_RP
         self % numTimeSteps = 0
         self % dt           = 0.0_RP
         
      END SUBROUTINE destructTimeIntegrator
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE Integrate( self, sem, controlVariables, monitors, pAdaptator, ComputeTimeDerivative, ComputeTimeDerivativeIsolated, &
                            ComputeTimeDerivative_onlyLinear, ComputeTimeDerivative_onlyNonLinear) 
      use pAdaptationClass
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
      type(pAdaptation_t)                  :: pAdaptator
      procedure(ComputeQDot_FCN)           :: ComputeTimeDerivative
      procedure(ComputeQDot_FCN)           :: ComputeTimeDerivativeIsolated
      procedure(ComputeQDot_FCN), optional :: ComputeTimeDerivative_onlyLinear
      procedure(ComputeQDot_FCN), optional :: ComputeTimeDerivative_onlyNonLinear

!
!     ---------
!     Internal variables
!     ---------
!
      integer              :: PA_Stage  ! P-adaptation stage
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
      
!     Perform FMG cycle if requested
!        (only for steady simulations)
!     ------------------------------
      
      if (self % integratorType == STEADY_STATE .and. &
          controlVariables % containsKey("fasfmg residual")) then
          
         FMGres = controlVariables % doubleprecisionValueForKey("fasfmg residual")
         write(STD_OUT,*) 'Using FMG solver to get initial condition. Res =', FMGres
         
         call FMGSolver % construct(controlVariables,sem)
         call FMGSolver % solve(0,0._RP, 0._RP, ComputeTimeDerivative, .TRUE.,FMGres) 
         
         call FMGSolver % destruct
      end if
      
!     Perform p-adaptation stage(s) if requested
!     ------------------------------------------
      if (pAdaptator % Adapt) then
         
         PA_Stage = 0
         do while (pAdaptator % Adapt)
            PA_Stage = PA_Stage + 1
            
            call IntegrateInTime( self, sem, controlVariables, monitors, ComputeTimeDerivative, pAdaptator % reqTE*0.1_RP)  ! The residual is hard-coded to 0.1 * truncation error threshold (see Kompenhans, Moritz, et al. "Adaptation strategies for high order discontinuous Galerkin methods based on Tau-estimation." Journal of Computational Physics 306 (2016): 216-236.)
            
            call pAdaptator % pAdaptTE(sem,sem  % numberOfTimeSteps,0._RP, ComputeTimeDerivative, ComputeTimeDerivativeIsolated, controlVariables)  ! Time is hardcoded to 0._RP (not important since it's only for STEADY_STATE)
            
            sem % numberOfTimeSteps = sem % numberOfTimeSteps + 1
            
         end do
      end if
      
!     Finish time integration
!     -----------------------
      if ( present(ComputeTimeDerivative_onlyLinear) ) then
         call IntegrateInTime( self, sem, controlVariables, monitors, ComputeTimeDerivative, CTD_linear = ComputeTimeDerivative_onlyLinear, CTD_nonlinear = ComputeTimeDerivative_onlyNonLinear)
      else
         call IntegrateInTime( self, sem, controlVariables, monitors, ComputeTimeDerivative)
      end if

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
   subroutine IntegrateInTime( self, sem, controlVariables, monitors, ComputeTimeDerivative, tolerance, CTD_linear, CTD_nonlinear)
      
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
      CLASS(TimeIntegrator_t)              :: self
      TYPE(DGSem)                          :: sem
      TYPE(FTValueDictionary), intent(in)  :: controlVariables
      class(Monitor_t)                     :: monitors
      procedure(ComputeQDot_FCN)           :: ComputeTimeDerivative
      real(kind=RP), optional, intent(in)  :: tolerance   !< ? tolerance to integrate down to
      procedure(ComputeQDot_FCN), optional :: CTD_linear
      procedure(ComputeQDot_FCN), optional :: CTD_nonlinear
!
!     ------------------
!     Internal variables
!     ------------------
!
interface
         subroutine UserDefinedPeriodicOperation(mesh, time, monitors)
            use SMConstants
            use HexMeshClass
            use MonitorsClass
            use PhysicsStorage
            IMPLICIT NONE
            CLASS(HexMesh)  :: mesh
            REAL(KIND=RP) :: time
            type(Monitor_t), intent(in)  :: monitors
         end subroutine UserDefinedPeriodicOperation
end interface
      
      real(kind=RP)                 :: Tol                                 ! Tolerance used for STEADY_STATE computations
      REAL(KIND=RP)                 :: t
      REAL(KIND=RP)                 :: maxResidual(NTOTALVARS)
      REAL(KIND=RP)                 :: dt
      INTEGER                       :: k, mNumber
      CHARACTER(LEN=13)             :: fName = "Movie_XX.tec"
      CHARACTER(LEN=2)              :: numChar
      CHARACTER(len=LINE_LENGTH)    :: SolutionFileName
      ! Time-step solvers:
      type(FASMultigrid_t)          :: FASSolver
      type(AnisFASMultigrid_t)      :: AnisFASSolver
      type(BDFIntegrator_t)         :: BDFSolver
      type(RosenbrockIntegrator_t)  :: RosenbrockSolver
      
      CHARACTER(len=LINE_LENGTH)    :: TimeIntegration
      logical                       :: saveGradients
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
      
      mNumber = 0
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
      call ComputeTimeDerivative(sem % mesh, sem % particles, t, sem % BCFunctions)
      maxResidual       = ComputeMaxResiduals(sem % mesh)
      sem % maxResidual = maxval(maxResidual)
      call Monitors % UpdateValues( sem % mesh, t, sem % numberOfTimeSteps, maxResidual )
      call self % Display(sem % mesh, monitors, sem  % numberOfTimeSteps)

      call monitors % WriteToFile(sem % mesh)

      IF (self % integratorType == STEADY_STATE) THEN
         IF (maxval(maxResidual) <= Tol )  THEN
            write(STD_OUT,'(/,A,I0,A,ES10.3)') "   *** Residual tolerance reached at iteration ",sem % numberOfTimeSteps," with Residual = ", maxval(maxResidual)
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
!        Autosave bounded time step
!        --------------------------
         dt = self % autosave % CorrectDt(t,self % dt)
!
!        Autosave bounded by time-accurate simulations
!        ---------------------------------------------
         if ( self % integratorType .eq. TIME_ACCURATE ) then
            if ( ( t + dt) .gt. self % tFinal ) then
               dt = self % tFinal - t
            end if
         end if
!
!        Perform time step
!        -----------------         
         SELECT CASE (TimeIntegration)
         CASE (IMPLICIT_SOLVER)
            call BDFSolver % TakeStep (sem, t , dt , ComputeTimeDerivative)
         CASE (ROSENBROCK_SOLVER)
            call RosenbrockSolver % TakeStep (sem, t , dt , ComputeTimeDerivative)
         CASE (EXPLICIT_SOLVER)
            CALL self % RKStep ( sem % mesh, sem % particles, t, sem % BCFunctions, dt, ComputeTimeDerivative)
         case (FAS_SOLVER)
            call FASSolver % solve(k, t, dt, ComputeTimeDerivative)
         case (ANISFAS_SOLVER)
            call AnisFASSolver % solve(k,t, ComputeTimeDerivative)
         case (IMEX_SOLVER)
            call TakeIMEXEulerStep(sem, t, dt, controlVariables, computeTimeDerivative, CTD_linear, CTD_nonlinear)
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
               write(STD_OUT,'(/,A,I0,A,ES10.3)') "   *** Residual tolerance reached at iteration ",k+1," with Residual = ", maxval(maxResidual)
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
!        User defined periodic operation
!        -------------------------------
         CALL UserDefinedPeriodicOperation(sem % mesh, t, monitors)
!
!        Print monitors
!        --------------
         IF( (MOD( k+1, self % outputInterval) == 0) .or. (k .eq. self % initial_iter) ) call self % Display(sem % mesh, monitors, k+1)
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
END MODULE TimeIntegratorClass
