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
      MODULE TimeIntegratorClass
      
      USE SMConstants
      use FTValueDictionaryClass
      USE PolynomialInterpAndDerivsModule
      USE DGSEMClass
      USE Physics
      USE ExplicitMethods
      use AutosaveClass
      use StopwatchClass
      IMPLICIT NONE 
      
      INTEGER, PARAMETER :: TIME_ACCURATE = 0, STEADY_STATE = 1

      TYPE TimeIntegrator_t
         INTEGER                                :: integratorType
         REAL(KIND=RP)                          :: tFinal, time, initial_time
         INTEGER                                :: initial_iter, numTimeSteps, outputInterval, iter
         REAL(KIND=RP)                          :: dt, tolerance, cfl
         LOGICAL                                :: Compute_dt                    ! Is st computed from an inputted CFL number?
         type(Autosave_t)                       :: autosave
         PROCEDURE(RKStepFcn), NOPASS , POINTER :: RKStep
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
            use DGSEMClass
            implicit none
            type(DGSem)     :: sem
            real(kind=RP)   :: t, deltaT
         end subroutine RKStepFcn
      end interface
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
            self % Compute_dt = .TRUE.
            self % cfl        = controlVariables % doublePrecisionValueForKey("cfl")
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
      SUBROUTINE Integrate( self, sem, controlVariables, monitors, pAdaptator)
      
      USE Implicit_JF , ONLY : TakeBDFStep_JF
      USE Implicit_NJ , ONLY : TakeBDFStep_NJ
      USE FASMultigridClass
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      CLASS(TimeIntegrator_t)       :: self
      TYPE(DGSem)                   :: sem
      TYPE(FTValueDictionary)       :: controlVariables
      class(Monitor_t)              :: monitors
      type(pAdaptation_t)           :: pAdaptator

!
!     ---------
!     Internal variables
!     ---------
!
      integer              :: PA_Stage  ! P-adaptation stage
      real(kind=RP)        :: FMGres    ! Target residual for FMG solver
      type(FASMultigrid_t) :: FMGSolver ! FAS multigrid solver for Full-Multigrid (FMG) initialization
      
!     Initializations
!     ---------------

      sem  % numberOfTimeSteps = self % initial_iter
      
!     Measure solver time
!     -------------------
      
      call Stopwatch % CreateNewEvent("Solver")
      call Stopwatch % Start("Solver")
      
!     Perform FMG cycle if requested
!     ------------------------------
      
      if (self % integratorType == STEADY_STATE .and. &
          controlVariables % containsKey("fasfmg residual")) then
          
         FMGres = controlVariables % RealValueForKey("fasfmg residual")
         write(STD_OUT,*) 'Using FMG solver to get initial condition. Res =', FMGres
         
         call FMGSolver % construct(controlVariables,sem)
         call FMGSolver % solve(0,0._RP,.TRUE.,FMGres)
         
         call FMGSolver % destruct
      end if
      
!     Perform p-adaptation stage(s) if requested
!     ------------------------------------------
      if (pAdaptator % Adapt) then
         write(STD_OUT,*) '*******    Performing p-adaptation    *******'
         
         PA_Stage = 0
         do while (pAdaptator % Adapt)
            PA_Stage = PA_Stage + 1
            
            call IntegrateInTime( self, sem, controlVariables, monitors, pAdaptator % reqTE*0.1_RP)  ! The residual is hard-coded to 0.1 * truncation error threshold (see Kompenhans, Moritz, et al. "Adaptation strategies for high order discontinuous Galerkin methods based on Tau-estimation." Journal of Computational Physics 306 (2016): 216-236.)
            
            !! TODO: Call p-Adaptator plotter
            call pAdaptator % pAdaptTE(sem,sem  % numberOfTimeSteps,0._RP)  ! Time is hardcoded to 0._RP (not important since it's only for STEADY_STATE)
            
            call self % Display(sem % mesh, monitors)
            
            !! TODO: Call p-Adaptator plotter
            !Write plot file
            
         end do
      end if
      
!     Finish time integration
!     -----------------------

      call IntegrateInTime( self, sem, controlVariables, monitors)

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
   subroutine IntegrateInTime( self, sem, controlVariables, monitors, tolerance)
      
      USE Implicit_JF , ONLY : TakeBDFStep_JF
      USE Implicit_NJ , ONLY : TakeBDFStep_NJ
      use FASMultigridClass
      use AnisFASMultigridClass
      use StopwatchClass
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      CLASS(TimeIntegrator_t)             :: self
      TYPE(DGSem)                         :: sem
      TYPE(FTValueDictionary), intent(in) :: controlVariables
      class(Monitor_t)                    :: monitors
      real(kind=RP), optional, intent(in) :: tolerance   !< ? tolerance to integrate down to

!
!     ------------------
!     Internal variables
!     ------------------
!
interface
         subroutine UserDefinedPeriodicOperation(mesh, time, monitors)
            use HexMeshClass
            use MonitorsClass
            IMPLICIT NONE
            CLASS(HexMesh)  :: mesh
            REAL(KIND=RP) :: time
            type(Monitor_t), intent(in)  :: monitors
         end subroutine UserDefinedPeriodicOperation
         character(len=LINE_LENGTH) function getFileName(inputLine)
            use SMConstants
            implicit none
            character(len=*)     :: inputLine
         end function getFileName
end interface
      
      real(kind=RP)                 :: Tol                                 ! Tolerance used for STEADY_STATE computations
      REAL(KIND=RP)                 :: t
      REAL(KIND=RP)                 :: maxResidual(N_EQN)
      REAL(KIND=RP)                 :: dt
      INTEGER                       :: k, mNumber
      CHARACTER(LEN=13)             :: fName = "Movie_XX.tec"
      CHARACTER(LEN=2)              :: numChar
      EXTERNAL                      :: ExternalState, ExternalGradients
      CHARACTER(len=LINE_LENGTH)    :: SolutionFileName
      ! For Implicit
      CHARACTER(len=LINE_LENGTH)    :: TimeIntegration
      INTEGER                       :: JacFlag
      type(FASMultigrid_t)          :: FASSolver
      type(AnisFASMultigrid_t)      :: AnisFASSolver
      logical                       :: saveGradients
!
!     ----------------------
!     Read Control variables
!     ----------------------
!
      IF (controlVariables % containsKey("time integration")) THEN
         TimeIntegration  = controlVariables % StringValueForKey("time integration",LINE_LENGTH)
      ELSE ! Default value
         TimeIntegration = 'explicit'
      END IF
      SolutionFileName   = trim(getFileName(controlVariables % StringValueForKey("solution file name",LINE_LENGTH)))
      
      ! Specific keywords
      IF (TimeIntegration == 'implicit') JacFlag = controlVariables % IntegerValueForKey("jacobian flag")
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
      
      if (TimeIntegration == 'FAS') CALL FASSolver % construct(controlVariables,sem)
      if (TimeIntegration == 'AnisFAS') CALL AnisFASSolver % construct(controlVariables,sem)
!
!     ------------------
!     Configure restarts
!     ------------------
!
      saveGradients = controlVariables % logicalValueForKey("save gradients with solution")
!
!     -----------------
!     Integrate in time
!     -----------------
!      
      mNumber = 0
      t = self % time
      sem % MaxResidual = 1.e-3_RP !initializing to this value for implicit solvers (Newton tolerance is computed according to this)
      
      DO k = self % initial_iter, self % initial_iter + self % numTimeSteps-1
!
!        CFL-bounded time step
!        ---------------------      
         IF ( self % Compute_dt ) self % dt = MaxTimeStep( sem, self % cfl )
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
            CASE ('implicit')
               SELECT CASE (JacFlag)
                  CASE (1)
                     CALL TakeBDFStep_JF (sem, t , dt )
                  CASE (2)
                     CALL TakeBDFStep_NJ (sem, t , dt , controlVariables)
                  CASE (3)
                     STOP 'Analytical Jacobian not implemented yet'
                  CASE DEFAULT
                     PRINT*, "Not valid 'Jacobian Flag'. Running with Jacobian-Free Newton-Krylov."
                     JacFlag = 1
                     CALL TakeBDFStep_JF (sem, t , dt )
               END SELECT
            CASE ('explicit')
               CALL self % RKStep ( sem, t, dt )
            case ('FAS')
               call FASSolver % solve(k,t)
            case ('AnisFAS')
               call AnisFASSolver % solve(k,t)
         END SELECT
!
!        Compute the new time
!        --------------------         
         t = t + dt
         self % time = t
!
!        Get maximum residuals
!        ---------------------
         maxResidual       = ComputeMaxResidual(sem)
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
              call self % Display(sem % mesh, monitors)
              sem  % maxResidual       = maxval(maxResidual)
              self % time              = t
              sem  % numberOfTimeSteps = k + 1

              write(STD_OUT,'(/,A,I0,A,ES10.3)') "   *** Residual tolerance reached at iteration ",k+1," with Residual = ", maxval(maxResidual)
              call Stopwatch % Pause("Solver")
              RETURN
            END IF
         ELSEIF (self % integratorType == TIME_ACCURATE) THEN
            IF ( t .ge. self % tFinal) then
               self % time = t     
               sem % numberOfTimeSteps = k+1
               call self % Display( sem % mesh, monitors)
               call Stopwatch % Pause("Solver")
               exit
            end if
         END IF
!
!        User defined periodic operation
!        -------------------------------
         CALL UserDefinedPeriodicOperation(sem % mesh, t, monitors)
!
!        Print monitors
!        --------------
         IF( (MOD( k+1, self % outputInterval) == 0) .or. (k .eq. self % initial_iter) ) call self % Display(sem % mesh, monitors)
!
!        Autosave
!        --------         
         if ( self % autosave % Autosave(k+1) ) then
            call SaveRestart(sem,k+1,t,SolutionFileName, saveGradients)
   
         end if

!        Flush monitors
!        --------------
         call monitors % WriteToFile(sem % mesh)

      END DO
!
!     Flush the remaining information in the monitors
!     -----------------------------------------------
      if ( k .ne. 0 ) then
         call Monitors % writeToFile(sem % mesh, force = .true. )
      end if
      
      sem % maxResidual       = maxval(maxResidual)
      self % time             = t
      sem % numberOfTimeSteps = k
      
!
!     ---------
!     Finish up
!     ---------
!
      if (TimeIntegration == 'FAS') CALL FASSolver % destruct
      if (TimeIntegration == 'AnisFAS') CALL AnisFASSolver % destruct
      
   end subroutine IntegrateInTime
      
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
!     Subroutine to print the residuals
!
!
   subroutine TimeIntegrator_Display(self, mesh, monitors)
      use MPI_Process_Info
      implicit none
      class(TimeIntegrator_t),   intent(in)     :: self
      class(HexMesh),            intent(in)     :: mesh
      class(Monitor_t),          intent(inout)  :: monitors
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
         if ( (self % integratorType .eq. TIME_ACCURATE) .and. (self % iter .gt. self % initial_iter+1) ) then 
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
      write(STD_OUT,'(A,A,A,ES10.3,A)') '*** Writing file "',trim(FinalName),'", with t = ',t,'.'
      call sem % mesh % SaveSolution(k,t,trim(finalName),saveGradients)
   
   END SUBROUTINE SaveRestart
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!         
END MODULE TimeIntegratorClass
