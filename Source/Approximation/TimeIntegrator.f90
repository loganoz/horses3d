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
      USE PolynomialInterpAndDerivsModule
      USE DGSEMClass
      USE Physics
      USE DGSEMPlotterClass
      USE UserDefinedFunctions
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
      END TYPE TimeIntegrator_t

      abstract interface
         subroutine RKStepFcn( sem , t , deltaT , maxResidual )
            use DGSEMClass
            implicit none
            type(DGSem)     :: sem
            real(kind=RP)   :: t, deltaT,  maxResidual(N_EQN)
         end subroutine RKStepFcn
      end interface
!
!     ========      
      CONTAINS 
!     ========      
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE constructTimeIntegrator(self,sem,controlVariables)
      
         IMPLICIT NONE
         !---------------------------------------------
         CLASS(TimeIntegrator_t)     :: self
         TYPE( DGSem )               :: sem
         TYPE(FTValueDictionary)     :: controlVariables
         !---------------------------------------------
         REAL(KIND=RP)               :: dt, cfl
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
         SELECT CASE (controlVariables % StringValueForKey("time integration",LINE_LENGTH))
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
      SUBROUTINE Integrate( self, sem, controlVariables)
      
      USE Implicit_JF , ONLY : TakeBDFStep_JF
      USE Implicit_NJ , ONLY : TakeBDFStep_NJ
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      CLASS(TimeIntegrator_t)       :: self
      TYPE(DGSem)                   :: sem
      TYPE(FTValueDictionary)       :: controlVariables

!
!     ---------
!     Internal variables
!     ---------
!
      
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
      LOGICAL                       :: imp !implicit?
      INTEGER                       :: JacFlag
!
!     ----------------------
!     Read Control variables
!     ----------------------
!
      imp              = controlVariables % LogicalValueForKey("implicit time")
      IF (imp) JacFlag = controlVariables % IntegerValueForKey("jacobian flag")
      RestFileName     = controlVariables % StringValueForKey("restart file name",LINE_LENGTH)
      RestartInterval  = controlVariables % IntegerValueForKey("restart interval") !If not present, RestartInterval=HUGE
!
!     -----------------
!     Integrate in time
!     -----------------
!      
      mNumber = 0
      t = self % time
      
      DO k = 0, self % numTimeSteps-1
      
         IF ( self % Compute_dt ) self % dt = MaxTimeStep( sem, self % cfl )
         
         IF (imp) THEN
            SELECT CASE (JacFlag)
               CASE (1)
                  CALL TakeBDFStep_JF (sem, t , self%dt , maxResidual)
               CASE (2)
                  CALL TakeBDFStep_NJ (sem, t , self%dt , maxResidual, controlVariables)
               CASE (3)
                  STOP 'Analytical Jacobian not implemented yet'
               CASE DEFAULT
                  PRINT*, "Not valid 'Jacobian Flag'. Running with Jacobian-Free Newton-Krylov."
                  JacFlag = 1
                  CALL TakeBDFStep_JF (sem, t , self%dt , maxResidual)
            END SELECT
         ELSE
            CALL self % RKStep ( sem, t, self % dt, maxResidual )
         END IF
         
         t = t + self % dt                    !arueda: changed since time step does not have to be constant!
         
         IF (self % integratorType == STEADY_STATE) THEN
            IF (maxval(maxResidual) <= self % tolerance )  THEN
              CALL PlotResiduals(k+1 , t , maxResidual)
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
               CALL PlotResiduals(k+1 , t , maxResidual)
               
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
        
      END DO
      
      sem % maxResidual       = maxval(maxResidual)
      self % time             = t
      sem % numberOfTimeSteps = k

      END SUBROUTINE Integrate
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE TakeRK3Step( sem, t, deltaT, maxResidual )
!
!     ----------------------------------
!     Williamson's 3rd order Runge-Kutta
!     ----------------------------------
!
      IMPLICIT NONE
!
!     -----------------
!     Input parameters:
!     -----------------
!
      TYPE(DGSem)     :: sem
      REAL(KIND=RP)   :: t, deltaT, tk, maxResidual(N_EQN)
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(3) :: a = (/0.0_RP       , -5.0_RP /9.0_RP , -153.0_RP/128.0_RP/)
      REAL(KIND=RP), DIMENSION(3) :: b = (/0.0_RP       ,  1.0_RP /3.0_RP ,    3.0_RP/4.0_RP/)
      REAL(KIND=RP), DIMENSION(3) :: c = (/1.0_RP/3.0_RP,  15.0_RP/16.0_RP,    8.0_RP/15.0_RP/)
      
      INTEGER :: k, id , eq
      REAL(KIND=RP) :: localMaxResidual(N_EQN)
!
      do id = 1, SIZE( sem % mesh % elements ) 
         sem % mesh % elements(id) % G = 0.0_RP   
      enddo 

      DO k = 1,3

         tk = t + b(k)*deltaT
         CALL ComputeTimeDerivative( sem, tk )

!$omp parallel do
         DO id = 1, SIZE( sem % mesh % elements )
            sem % mesh % elements(id) % G = a(k)*sem % mesh % elements(id) % G  +             sem % mesh % elements(id) % QDot
            sem % mesh % elements(id) % Q =      sem % mesh % elements(id) % Q  + c(k)*deltaT*sem % mesh % elements(id) % G
         END DO
!$omp end parallel do

      END DO
!
!     ----------------
!     Compute residual
!     ----------------
!
      maxResidual = 0.0_RP
      DO id = 1, SIZE( sem % mesh % elements )
         DO eq = 1 , N_EQN
            localMaxResidual(eq) = MAXVAL(ABS(sem % mesh % elements(id) % QDot(:,:,:,eq)))
            maxResidual(eq) = MAX(maxResidual(eq),localMaxResidual(eq))
         END DO
      END DO
      
   END SUBROUTINE TakeRK3Step
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
!     Subroutine to print the residuals
!
!
   subroutine PlotResiduals( iter , time , maxResiduals )
      implicit none
      integer, intent(in)       :: iter
      real(kind=RP), intent(in) :: time
      real(kind=RP), intent(in) :: maxResiduals(N_EQN)
!     --------------------------------------------------------
      integer, parameter        :: showLabels = 50
      integer, save             :: shown = 0

      if ( mod(shown , showLabels) .eq. 0 ) then ! Show labels
         write(STD_OUT , '(/)')
         write(STD_OUT , '(/)')
         write(STD_OUT , '(A17,3X,A10,3X,A10,3X,A10,3X,A10,3X,A10,3X,A10)') &
               "Iteration" , "time" , "continuity" , "x-momentum" , "y-momentum", &
               "z-momentum" , "energy"
         write(STD_OUT , '(A17,3X,A10,3X,A10,3X,A10,3X,A10,3X,A10,3X,A10)') &
               "---------" , "--------" , "----------" , "----------" , "----------" , &
               "----------", "--------"
      end if
      write(STD_OUT , 110) iter ,"|", time ,"|", maxResiduals(IRHO) , "|" , maxResiduals(IRHOU) , &
                                          "|", maxResiduals(IRHOV) , "|" , maxResiduals(IRHOW) , "|" , maxResiduals(IRHOE)
      110 format (I17,X,A,X,ES10.3,X,A,X,ES10.3,X,A,X,ES10.3,X,A,X,ES10.3,X,A,X,ES10.3,X,A,X,ES10.3)

    shown = shown + 1

   end subroutine PlotResiduals
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
