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
      USE PDEModule
      IMPLICIT NONE 
!
      INTEGER, PARAMETER :: TIME_ACCURATE = 0, STEADY_STATE = 1
      
      TYPE TimeIntegrator
         INTEGER       :: integratorType
         REAL(KIND=RP) :: tFinal, tStart
         INTEGER       :: numTimeSteps, plotInterval
         REAL(KIND=RP) :: dt, tolerance
      END TYPE TimeIntegrator
!
!     --------
!     Generics
!     --------
!
      INTERFACE NewTimeIntegrator
         MODULE PROCEDURE NewAccurateTimeIntegrator
         MODULE PROCEDURE NewSteadyTimeIntegrator
      END INTERFACE NewTimeIntegrator
      INTERFACE Destruct
         MODULE PROCEDURE DestructTimeIntegrator
      END INTERFACE Destruct
!
!     ========      
      CONTAINS 
!     ========      
!
      TYPE(TimeIntegrator) FUNCTION NewAccurateTimeIntegrator( startTime, finalTime, numberOfSteps, plotInterval )  RESULT(this)
         IMPLICIT NONE
         REAL(KIND=RP)              :: finalTime, startTime
         INTEGER                    :: numberOfSteps, plotInterval
!
!        --------------------------------------------------------------
!        Compute time step and set values for time accurate computation
!        --------------------------------------------------------------
!
         this % tFinal         = finalTime
         this % tStart         = startTime
         this % numTimeSteps   = numberOfSteps
         this % dt             = (this % tFinal - this % tStart)/numberOfSteps
         this % plotInterval   = plotInterval
         this % integratorType = TIME_ACCURATE
         this % tolerance      = 1.d-11
      
      END FUNCTION NewAccurateTimeIntegrator
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      TYPE(TimeIntegrator) FUNCTION NewSteadyTimeIntegrator( dt, numberOfSteps, plotInterval )  RESULT(this)
         IMPLICIT NONE
         REAL(KIND=RP)              :: dt
         INTEGER                    :: numberOfSteps, plotInterval
!
!        ---------------------------------------------------------------
!        Compute time step and set values for a steady-state computation
!        ---------------------------------------------------------------
!
         this % numTimeSteps   = numberOfSteps
         this % dt             = dt
         this % plotInterval   = plotInterval
         this % integratorType = STEADY_STATE
         this % tolerance      = 1.d-11
      
      END FUNCTION NewSteadyTimeIntegrator
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructTimeIntegrator( this ) 
         TYPE(TimeIntegrator) :: this
         this % tFinal       = 0.0_RP
         this % tStart       = 0.0_RP
         this % numTimeSteps = 0
         this % dt           = 0.0_RP
      END SUBROUTINE DestructTimeIntegrator
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SetIterationTolerance( this, tol ) 
         TYPE(TimeIntegrator) :: this
         REAL(KIND=RP)        :: tol
         this % tolerance = tol
      END SUBROUTINE SetIterationTolerance
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      REAL(KIND=RP) FUNCTION MaxTimeStep( sem, cfl ) 
         IMPLICIT NONE
         TYPE(DGSem)    :: sem
         REAL(KIND=RP)  :: cfl
         REAL(KIND=RP), EXTERNAL :: MaximumEigenvalue ! Supplied in the PDEModule
         
         MaxTimeStep  = cfl/MaximumEigenvalue( sem )
      
      END FUNCTION MaxTimeStep
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE Integrate( this, sem, cfl, thePlotter, ExternalState, ExternalGradients )
      
      USE PlotterClass
      IMPLICIT NONE
      
      TYPE(TimeIntegrator)  :: this
      TYPE(DGSem)           :: sem
      TYPE(Plotter)         :: thePlotter
      REAL(KIND=RP)         :: cfl
      
      REAL(KIND=RP)         :: t, maxResidual
      INTEGER               :: k, pUnit = 11, mNumber
      CHARACTER(LEN=13)     :: fName = "Movie_XX.tec"
      CHARACTER(LEN=2)      :: numChar
      EXTERNAL              :: ExternalState, ExternalGradients
!
!     -----------------
!     Integrate in time
!     -----------------
!
      mNumber = 0
      DO k = 0, this % numTimeSteps-1
         IF ( this % integratorType == STEADY_STATE ) THEN
            this % dt = MaxTimeStep( sem, cfl )
         END IF

         t = this % tStart + k*this % dt
         CALL TakeRK3Step( sem, t, this % dt, ExternalState, ExternalGradients, maxResidual )
         IF( this % integratorType == STEADY_STATE .AND. maxResidual <= this % tolerance )     THEN
            PRINT *, "Residual tolerance reached. Residual = ", maxResidual
            RETURN
         END IF

         IF( MOD( k+1, this % plotInterval) == 0 )     THEN
            IF ( this % integratorType == STEADY_STATE )     THEN
               PRINT *, k, CHAR(9), LOG10(maxResidual)
            ELSE
               mNumber = mNumber + 1
               WRITE(numChar,'(i2)') mNumber
               IF ( mNumber >= 10 )     THEN
                  fName(7:8) = numChar
               ELSE
                  fName(7:7) = "0"
                  fName(8:8) = numChar(2:2)
               END IF
               PRINT *, fName
               OPEN(UNIT=thePlotter % fUnit, FILE = fName)                           
                  CALL ExportToTecplot( thePlotter, sem )               
               CLOSE(thePlotter % fUnit)
            END IF
        END IF
      END DO

      END SUBROUTINE Integrate
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE TakeRK3Step( sem, t, deltaT, ExternalState, ExternalGradients, maxResidual )
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
      
      INTEGER :: k, id, i, j, m
      REAL(KIND=RP) :: localMaxResidual
!
      
      DO k = 1,3
         tk = t + b(k)*deltaT
         CALL ComputeTimeDerivative( sem, tk, ExternalState, ExternalGradients )
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
!$omp parallel do
         DO id = 1, SIZE( sem % mesh % elements )
            sem % mesh % elements(id) % G = a(k)*sem % mesh % elements(id) % G  +             sem % mesh % elements(id) % QDot
            sem % mesh % elements(id) % Q =       sem % mesh % elements(id) % Q  + c(k)*deltaT*sem % mesh % elements(id) % G
         END DO
!$omp end parallel do
      END DO
      
   END SUBROUTINE TakeRK3Step

END MODULE TimeIntegratorClass
