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
      USE DGSEMPlotterClass
      USE ExplicitMethods
      IMPLICIT NONE 
!
      INTEGER, PARAMETER :: TIME_ACCURATE = 0, STEADY_STATE = 1
      
      TYPE TimeIntegrator_t
         INTEGER                                :: integratorType
         REAL(KIND=RP)                          :: tFinal, time
         INTEGER                                :: numTimeSteps, plotInterval
         REAL(KIND=RP)                          :: dt, tolerance, cfl
         LOGICAL                                :: Compute_dt                    ! Is st computed from an inputted CFL number?
         TYPE(DGSEMPlotter),   POINTER          :: plotter  !Plotter is NOT owned by the time integrator
         PROCEDURE(RKStepFcn), NOPASS , POINTER :: RKStep
!
!        ========         
         CONTAINS
!        ========         
!
         PROCEDURE :: construct => constructTimeIntegrator
         PROCEDURE :: destruct => destructTimeIntegrator
         PROCEDURE :: setPlotter
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
      SUBROUTINE constructTimeIntegrator(self,controlVariables)
      
         IMPLICIT NONE
         !---------------------------------------------
         CLASS(TimeIntegrator_t)     :: self
         TYPE(FTValueDictionary)     :: controlVariables
         !---------------------------------------------
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
         self % time           =  0._RP                                                                  ! TODO: Modify this for restarted cases?
         self % numTimeSteps   =  controlVariables % integerValueForKey ("number of time steps")
         self % plotInterval   =  controlVariables % integerValueForKey("output interval")
         self % tolerance      =  controlVariables % doublePrecisionValueForKey("convergence tolerance")
         self % plotter        => NULL()
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
      SUBROUTINE setPlotter( self, plotter ) 
         CLASS(TimeIntegrator_t)     :: self
         TYPE(DGSEMPlotter), pointer :: plotter
         self % plotter => plotter
      END SUBROUTINE setPlotter
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE Integrate( self, sem, controlVariables, monitors)
      
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

!
!     ---------
!     Internal variables
!     ---------
!
interface
         SUBROUTINE UserDefinedPeriodicOperation(sem, time)
            USE DGSEMClass
            IMPLICIT NONE
            CLASS(DGSem)  :: sem
            REAL(KIND=RP) :: time
         END SUBROUTINE UserDefinedPeriodicOperation
end interface
      
      REAL(KIND=RP)                 :: t
      REAL(KIND=RP)                 :: maxResidual(N_EQN)
      INTEGER                       :: k, mNumber
      CHARACTER(LEN=13)             :: fName = "Movie_XX.tec"
      CHARACTER(LEN=2)              :: numChar
      EXTERNAL                      :: ExternalState, ExternalGradients
      
      ! For saving restarts
      CHARACTER(len=LINE_LENGTH)    :: RestFileName
      INTEGER                       :: RestartInterval
      
      ! For Implicit
      CHARACTER(len=LINE_LENGTH)    :: TimeIntegration
      INTEGER                       :: JacFlag
      
      TYPE(FASMultigrid_t)          :: FASSolver
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
      RestFileName     = controlVariables % StringValueForKey("restart file name",LINE_LENGTH)
      RestartInterval  = controlVariables % IntegerValueForKey("restart interval") !If not present, RestartInterval=HUGE
      
      ! Specific keywords
      IF (TimeIntegration == 'implicit') JacFlag = controlVariables % IntegerValueForKey("jacobian flag")
!
!     ---------------
!     Initializations
!     ---------------
!
      IF (TimeIntegration == 'FAS') CALL FASSolver % construct(controlVariables,sem)
!
!     -----------------
!     Integrate in time
!     -----------------
!      
      mNumber = 0
      t = self % time
      sem % MaxResidual = 1.e-3_RP !initializing to this value for implicit solvers (Newton tolerance is computed according to this)
      
      DO k = 0, self % numTimeSteps-1
      
         IF ( self % Compute_dt ) self % dt = MaxTimeStep( sem, self % cfl )
!
!        Perform time step
!        -----------------         
         SELECT CASE (TimeIntegration)
            CASE ('implicit')
               SELECT CASE (JacFlag)
                  CASE (1)
                     CALL TakeBDFStep_JF (sem, t , self%dt )
                  CASE (2)
                     CALL TakeBDFStep_NJ (sem, t , self%dt , controlVariables)
                  CASE (3)
                     STOP 'Analytical Jacobian not implemented yet'
                  CASE DEFAULT
                     PRINT*, "Not valid 'Jacobian Flag'. Running with Jacobian-Free Newton-Krylov."
                     JacFlag = 1
                     CALL TakeBDFStep_JF (sem, t , self%dt )
               END SELECT
            CASE ('explicit')
               CALL self % RKStep ( sem, t, self % dt )
            CASE ('FAS')
               CALL FASSolver % solve(k,t)
         END SELECT
!
!        Compute the new time
!        --------------------         
         t = t + self % dt
!
!        Get maximum residuals
!        ---------------------
         maxResidual       = ComputeMaxResidual(sem)
         sem % maxResidual = maxval(maxResidual)
!
!        Exit if the time exceeds the final time (only in TIME_ACCURATE mode)
!        ---------------------------------------
         if ( self % integratorType .eq. TIME_ACCURATE ) then
            if ( self % time .ge. self % tFinal) then
               call monitors % UpdateValues( sem % mesh, t, k+1, maxResidual )
               call self % Display( sem % mesh, monitors)
               exit
            end if
         end if
!
!        Update monitors
!        ---------------
         call Monitors % UpdateValues( sem % mesh, t, k+1, maxResidual )

         IF (self % integratorType == STEADY_STATE) THEN
            IF (maxval(maxResidual) <= self % tolerance )  THEN
              call self % Display(sem % mesh, monitors)
              sem % maxResidual       = maxval(maxResidual)
              self % time             = t
              sem % numberOfTimeSteps = k + 1
              PRINT *, "Residual tolerance reached at iteration ",k+1," with Residual = ", maxResidual
              RETURN
            END IF
         ELSEIF (self % integratorType == TIME_ACCURATE) THEN
            self % time              = t     
            IF (self % time > self % tFinal) EXIT
         END IF
         
         IF (MOD( k+1, RestartInterval) == 0) CALL SaveRestart(sem,k+1,t,RestFileName)
         
         IF( (MOD( k+1, self % plotInterval) == 0) .or. (k .eq. 0) )     THEN
          CALL UserDefinedPeriodicOperation(sem,t)
          
            IF ( self % integratorType == STEADY_STATE )     THEN
               call self % Display(sem % mesh, monitors)
               
            ELSE IF (ASSOCIATED(self % plotter))     THEN 
               mNumber = mNumber + 1
               WRITE(numChar,'(i2)') mNumber
               IF ( mNumber >= 10 )     THEN
                  fName(7:8) = numChar
               ELSE
                  fName(7:7) = "0"
                  fName(8:8) = numChar(2:2)
               END IF
               PRINT *, fName
               OPEN(UNIT = self % plotter % fUnit, FILE = fName)
                  CALL self % plotter % ExportToTecplot( sem % mesh % elements, sem % spA )
               CLOSE(self % plotter % fUnit)
               
            END IF
         END IF
!
!        Flush monitors
!        --------------
         call monitors % WriteToFile()
        
      END DO

      if ( k .ne. 0 ) then
         call Monitors % writeToFile( force = .true. )
      end if
      
      sem % maxResidual       = maxval(maxResidual)
      self % time             = t
      sem % numberOfTimeSteps = k

      END SUBROUTINE Integrate    
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
!     Subroutine to print the residuals
!
!
   subroutine TimeIntegrator_Display(self, mesh, monitors)
      implicit none
      class(TimeIntegrator_t),   intent(in)     :: self
      class(HexMesh),            intent(in)     :: mesh
      class(Monitor_t),          intent(inout)  :: monitors
!
!     ---------------
!     Local variables      
!     ---------------
!
      integer, parameter      :: showLabels = 50
      integer, save           :: shown = 0

      if ( mod(shown, showLabels) .eq. 0 ) then
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
   SUBROUTINE SaveRestart(sem,k,t,RestFileName)
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
!     ----------------------------------------------
      INTEGER                      :: fd             !  File unit for new restart file
      CHARACTER(len=LINE_LENGTH)   :: FinalName      !  Final name for particular restart file
!     ----------------------------------------------
      
      WRITE(FinalName,'(2A,I5.5,A,ES9.3,A)')  TRIM(RestFileName),'_step_',k,'_time_',t,'.rst'
      
      OPEN( newunit = fd             , &
            FILE    = TRIM(FinalName), & 
            ACTION  = 'WRITE'        , &
            FORM = "UNFORMATTED")
         CALL SaveSolutionForRestart( sem, fd )
      CLOSE (fd)
   
   END SUBROUTINE SaveRestart
END MODULE TimeIntegratorClass
