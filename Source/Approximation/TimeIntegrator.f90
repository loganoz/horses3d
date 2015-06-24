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
      IMPLICIT NONE 
!
      INTEGER, PARAMETER :: TIME_ACCURATE = 0, STEADY_STATE = 1
      
      TYPE TimeIntegrator
         INTEGER                     :: integratorType
         REAL(KIND=RP)               :: tFinal, tStart
         INTEGER                     :: numTimeSteps, plotInterval
         REAL(KIND=RP)               :: dt, tolerance, cfl
         TYPE(DGSEMPlotter), POINTER :: plotter
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
      END TYPE TimeIntegrator
!
!     ========      
      CONTAINS 
!     ========      
!
      SUBROUTINE constructAsTimeAccurateIntegrator( self,  startTime, finalTime, numberOfSteps, plotInterval )
         IMPLICIT NONE
         CLASS(TimeIntegrator) :: self
         REAL(KIND=RP)         :: finalTime, startTime
         INTEGER               :: numberOfSteps, plotInterval
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
      
      END SUBROUTINE constructAsTimeAccurateIntegrator
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE constructAsSteadyStateIntegrator( self, dt, cfl, numberOfSteps, plotInterval )
         IMPLICIT NONE
         CLASS(TimeIntegrator) :: self
         REAL(KIND=RP)         :: dt
         REAL(KIND=RP)         :: cfl
         INTEGER               :: numberOfSteps, plotInterval
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
      
      END SUBROUTINE constructAsSteadyStateIntegrator
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE destructTimeIntegrator( self ) 
         CLASS(TimeIntegrator) :: self
         self % tFinal       = 0.0_RP
         self % tStart       = 0.0_RP
         self % numTimeSteps = 0
         self % dt           = 0.0_RP
         
         IF(ASSOCIATED(self % plotter))     THEN
            CALL self % plotter % Destruct()
            DEALLOCATE(self % plotter)
         END IF 
         
      END SUBROUTINE destructTimeIntegrator
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE setPlotter( self, plotter ) 
         CLASS(TimeIntegrator)       :: self
         TYPE(DGSEMPlotter), pointer :: plotter
         self % plotter => plotter
      END SUBROUTINE setPlotter
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE setIterationTolerance( self, tol ) 
         CLASS(TimeIntegrator) :: self
         REAL(KIND=RP)         :: tol
         self % tolerance = tol
      END SUBROUTINE setIterationTolerance
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      REAL(KIND=RP) FUNCTION MaxTimeStep( sem, cfl ) 
         IMPLICIT NONE
         TYPE(DGSem)    :: sem
         REAL(KIND=RP)  :: cfl
         
         MaxTimeStep  = cfl/MaximumEigenvalue( sem )
      
      END FUNCTION MaxTimeStep
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
      TYPE(TimeIntegrator)  :: self
      TYPE(DGSem)           :: sem
      
      REAL(KIND=RP)         :: t, maxResidual
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
         
         CALL TakeRK3Step( sem, t, self % dt, maxResidual )
         
         IF( self % integratorType == STEADY_STATE .AND. maxResidual <= self % tolerance )     THEN
            PRINT *, "Residual tolerance reached. Residual = ", maxResidual
            RETURN
         END IF
         
         IF( MOD( k+1, self % plotInterval) == 0 )     THEN
         
            IF ( self % integratorType == STEADY_STATE )     THEN
               PRINT *, k, CHAR(9), LOG10(maxResidual)
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
      REAL(KIND=RP)   :: t, deltaT, tk, maxResidual
      EXTERNAL        :: ExternalState, ExternalGradients
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(3) :: a = (/0.0_RP       , -5.0_RP /9.0_RP , -153.0_RP/128.0_RP/)
      REAL(KIND=RP), DIMENSION(3) :: b = (/0.0_RP       ,  1.0_RP /3.0_RP ,    3.0_RP/4.0_RP/)
      REAL(KIND=RP), DIMENSION(3) :: c = (/1.0_RP/3.0_RP,  15.0_RP/16.0_RP,    8.0_RP/15.0_RP/)
      
      INTEGER :: k, id
      REAL(KIND=RP) :: localMaxResidual
!
      
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
               localMaxResidual = MAXVAL(ABS(sem % mesh % elements(id) % QDot))
               maxResidual = MAX(maxResidual,localMaxResidual)
            END DO
         END IF

         DO id = 1, SIZE( sem % mesh % elements )
            sem % mesh % elements(id) % G = a(k)*sem % mesh % elements(id) % G  +             sem % mesh % elements(id) % QDot
            sem % mesh % elements(id) % Q =      sem % mesh % elements(id) % Q  + c(k)*deltaT*sem % mesh % elements(id) % G
         END DO

      END DO
      
   END SUBROUTINE TakeRK3Step

END MODULE TimeIntegratorClass
