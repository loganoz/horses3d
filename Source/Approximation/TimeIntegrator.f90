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
      
      TYPE RKTimeIntegrator
         INTEGER                                :: integratorType
         REAL(KIND=RP)                          :: tFinal, tStart, time
         INTEGER                                :: numTimeSteps, plotInterval
         REAL(KIND=RP)                          :: dt, tolerance, cfl
         TYPE(DGSEMPlotter),   POINTER          :: plotter  !Plotter is NOT owned by the time integrator
         PROCEDURE(RKStepFcn), NOPASS , POINTER :: RKStep
!
!        ========         
         CONTAINS
!        ========         
!
         PROCEDURE :: constructAsTimeAccurateIntegrator
         PROCEDURE :: constructAsSteadyStateIntegrator
         PROCEDURE :: destruct => destructTimeIntegrator
         PROCEDURE :: setIterationTolerance
         PROCEDURE :: setPlotter
         PROCEDURE :: integrate
      END TYPE RKTimeIntegrator

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
      SUBROUTINE constructAsTimeAccurateIntegrator( self,  startTime, finalTime, numberOfSteps, plotInterval )
         IMPLICIT NONE
         CLASS(RKTimeIntegrator) :: self
         REAL(KIND=RP)           :: finalTime, startTime
         INTEGER                 :: numberOfSteps, plotInterval
!
!        --------------------------------------------------------------
!        Compute time step and set values for time accurate computation
!        --------------------------------------------------------------
!
         self % tFinal         = finalTime
         self % tStart         = startTime
         self % numTimeSteps   = numberOfSteps
         self % dt             = (self % tFinal - self % tStart)/numberOfSteps
         self % plotInterval   = plotInterval
         self % integratorType = TIME_ACCURATE
         self % tolerance      = 1.d-11
         self % plotter        => NULL()
         self % RKStep         => TakeRK3Step
      
      END SUBROUTINE constructAsTimeAccurateIntegrator
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE constructAsSteadyStateIntegrator( self, dt, cfl, numberOfSteps, plotInterval )
         IMPLICIT NONE
         CLASS(RKTimeIntegrator) :: self
         REAL(KIND=RP)           :: dt
         REAL(KIND=RP)           :: cfl
         INTEGER                 :: numberOfSteps, plotInterval
!
!        ---------------------------------------------------------------
!        Compute time step and set values for a steady-state computation
!        ---------------------------------------------------------------
!
         self % numTimeSteps   = numberOfSteps
         self % dt             = dt
         self % plotInterval   = plotInterval
         self % integratorType = STEADY_STATE
         self % tolerance      = 1.d-11
         self % cfl            = cfl
         self % plotter        => NULL()
         self % RKStep         => TakeRK3Step
      
      END SUBROUTINE constructAsSteadyStateIntegrator
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE destructTimeIntegrator( self ) 
         CLASS(RKTimeIntegrator) :: self
         self % tFinal       = 0.0_RP
         self % tStart       = 0.0_RP
         self % numTimeSteps = 0
         self % dt           = 0.0_RP
         
      END SUBROUTINE destructTimeIntegrator
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE setPlotter( self, plotter ) 
         CLASS(RKTimeIntegrator)     :: self
         TYPE(DGSEMPlotter), pointer :: plotter
         self % plotter => plotter
      END SUBROUTINE setPlotter
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE setIterationTolerance( self, tol ) 
         CLASS(RKTimeIntegrator) :: self
         REAL(KIND=RP)           :: tol
         self % tolerance = tol
      END SUBROUTINE setIterationTolerance
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE Integrate( self, sem )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      CLASS(RKTimeIntegrator)  :: self
      TYPE(DGSem)             :: sem
      
      REAL(KIND=RP)         :: t
      REAL(KIND=RP)         :: maxResidual(N_EQN)
      INTEGER               :: k, mNumber
      CHARACTER(LEN=13)     :: fName = "Movie_XX.tec"
      CHARACTER(LEN=2)      :: numChar
      EXTERNAL              :: ExternalState, ExternalGradients
!
!     -----------------
!     Integrate in time
!     -----------------
!
      mNumber = 0
      
      DO k = 0, self % numTimeSteps-1
      
         IF ( self % integratorType == STEADY_STATE ) THEN
            self % dt = MaxTimeStep( sem, self % cfl )
         END IF

         t = self % tStart + k*self % dt
         
         CALL self % RKStep ( sem, t, self % dt, maxResidual )

         IF( self % integratorType == STEADY_STATE .AND. maxval(maxResidual) <= self % tolerance )     THEN
         
            sem % maxResidual       = maxval(maxResidual)
            self % time             = t
            sem % numberOfTimeSteps = k + 1
            PRINT *, "Residual tolerance reached at iteration ",k+1," with Residual = ", maxResidual
            RETURN
            
         END IF
         
         IF( (MOD( k+1, self % plotInterval) == 0) .or. (k .eq. 0) )     THEN

          CALL UserDefinedPeriodicOperation(sem,t)
            
            IF ( self % integratorType == STEADY_STATE )     THEN
               CALL PlotResiduals(k+1 , t+self % dt , maxResidual)
               
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
                  CALL self % plotter % ExportToTecplot( sem % mesh % elements )
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
!
!        ----------------
!        Compute residual
!        ----------------
!
         IF ( k == 1 )     THEN
            maxResidual = 0.0_RP
            DO id = 1, SIZE( sem % mesh % elements )
               DO eq = 1 , N_EQN
                  localMaxResidual(eq) = MAXVAL(ABS(sem % mesh % elements(id) % QDot(:,:,:,eq)))
                  maxResidual(eq) = MAX(maxResidual(eq),localMaxResidual(eq))
               END DO
            END DO
         END IF

!$omp parallel do
         DO id = 1, SIZE( sem % mesh % elements )
            sem % mesh % elements(id) % G = a(k)*sem % mesh % elements(id) % G  +             sem % mesh % elements(id) % QDot
            sem % mesh % elements(id) % Q =      sem % mesh % elements(id) % Q  + c(k)*deltaT*sem % mesh % elements(id) % G
         END DO
!$omp end parallel do

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
         write(STD_OUT , '(A20,5X,A10,5X,A10,5X,A10,5X,A10,5X,A10,5X,A10)') &
               "Iteration" , "time" , "continuity" , "x-momentum" , "y-momentum", &
               "z-momentum" , "energy"
         write(STD_OUT , '(A20,5X,A10,5X,A10,5X,A10,5X,A10,5X,A10,5X,A10)') &
               "---------" , "--------" , "----------" , "----------" , "----------" , &
               "----------", "--------"
      end if
      write(STD_OUT , 110) iter ,"|", time ,"|", maxResiduals(IRHO) , "|" , maxResiduals(IRHOU) , &
                                          "|", maxResiduals(IRHOV) , "|" , maxResiduals(IRHOW) , "|" , maxResiduals(IRHOE)
      110 format (I20,2X,A,2X,ES10.3,2X,A,2X,ES10.3,2X,A,2X,ES10.3,2X,A,2X,ES10.3,2X,A,2X,ES10.3,2X,A,2X,ES10.3)

    shown = shown + 1

   end subroutine PlotResiduals

END MODULE TimeIntegratorClass
