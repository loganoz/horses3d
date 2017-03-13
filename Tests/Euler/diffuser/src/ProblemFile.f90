!
!////////////////////////////////////////////////////////////////////////
!
!      ProblemFile.f90
!      Created: June 26, 2015 at 8:47 AM 
!      By: David Kopriva  
!
!      The Problem File contains user defined procedures
!      that are used to "personalize" i.e. define a specific
!      problem to be solved. These procedures include initial conditions,
!      exact solutions (e.g. for tests), etc. and allow modifications 
!      without having to modify the main code.
!
!      The procedures, *even if empty* that must be defined are
!
!      UserDefinedStartup
!      UserDefinedInitialCondition(sem)
!      UserDefinedPeriodicOperation(sem)
!      UserDefinedFinalize(sem)
!      UserDefinedTermination
!
!      *** This problem file sets up a subsonic point source *** 
!
!//////////////////////////////////////////////////////////////////////// 
!
      MODULE UserDefinedDataStorage
         USE SMConstants
         IMPLICIT NONE 
         REAL(KIND=RP) :: rad0, f, h 
      END MODULE UserDefinedDataStorage
!
!//////////////////////////////////////////////////////////////////////// 
! 
      MODULE UserDefinedFunctions

!
!     ========      
      CONTAINS
!     ========
!
         SUBROUTINE UserDefinedStartup  
!
!        --------------------------------
!        Called before any other routines
!        --------------------------------
!
            IMPLICIT NONE  
         END SUBROUTINE UserDefinedStartup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalSetup(sem, controlVariables)
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in but before time integration
!           to allow mesh related initializations or memory allocations
!           ----------------------------------------------------------------------
!
            USE DGSEMClass
            USE Physics
            USE UserDefinedDataStorage
            USE FTValueDictionaryClass
            
            IMPLICIT NONE
!
!           ---------
!           Arguments
!           ---------
!
            CLASS(DGSem)            :: sem
            TYPE(FTValueDictionary) :: controlVariables
!
!           ---------------
!           Local variables
!           ---------------
!
            INTEGER       :: nodeID
            REAL(KIND=RP) :: x(3)
            REAL(KIND=RP) :: rad
!
!           --------------------------------
!           Set up for the diffuser geometry
!           --------------------------------
!
            rad0 = HUGE(1.0_RP)
            DO nodeID = 1, SIZE(sem % mesh % nodes)
               x   = sem % mesh % nodes(nodeID) % x
               rad = SQRT(x(1)**2 + x(2)**2)
               rad0  = MIN(rad0, rad)
            END DO
            
            f = sqrtGamma*rad0*mach
            h = GAMMA *(1.0_RP/gammaMinus1Div2 + mach**2)
            
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedInitialCondition(sem, controlVariables)
!
!           ------------------------------------------------
!           Called to set the initial condition for the flow
!           ------------------------------------------------
!
            USE SMConstants
            USE DGSEMClass
            USE PhysicsStorage
            IMPLICIT NONE
            
            TYPE(DGSem)              :: sem
            class(FTValueDictionary) :: controlVariables
            LOGICAL                  :: success
                     
            INTEGER     :: i, j, k, eID
            
            DO eID = 1, SIZE(sem % mesh % elements)
               DO k = 0, sem % spA % N
                  DO j = 0, sem % spA % N
                     DO i = 0, sem % spA % N 
                        CALL pointSourceFlowSolution( sem % mesh % elements(eID) % geom % x(:,i,j,k), &
                                                      sem % mesh % elements(eID) % Q(i,j,k,1:N_EQN), success )
                        IF(.NOT. success) ERROR STOP "Unable to compute initial condition"       
                     END DO
                  END DO
               END DO 
            END DO 
            
         END SUBROUTINE UserDefinedInitialCondition
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedPeriodicOperation(sem, time)
!
!           ----------------------------------------------------------
!           Called at the output interval to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            USE DGSEMClass
            IMPLICIT NONE
            CLASS(DGSem)  :: sem
            REAL(KIND=RP) :: time
            
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(sem, time)
            USE FTAssertions
!
!           --------------------------------------------------------
!           Called after the solution computed to allow, for example
!           error tests to be performed
!           --------------------------------------------------------
!
            USE DGSEMClass
            IMPLICIT NONE
!
!           ---------
!           Arguments
!           ---------
!
            CLASS(DGSem)  :: sem
            REAL(KIND=RP) :: time
!
!           ---------------
!           Local variables
!           ---------------
!
            CHARACTER(LEN=29)                  :: testName           = "Diffuser flow tests"
            REAL(KIND=RP)                      :: maxError
            REAL(KIND=RP), ALLOCATABLE         :: QExpected(:,:,:,:)
            INTEGER                            :: eID
            INTEGER                            :: i, j, k, N
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
!
!           -----------------------------------------------------------------------------------------
!           Expected solutions. Inflow/Outflow on all boundaries 
!           -----------------------------------------------------------------------------------------
!
!
!           ------------------------------------------------
!           Expected Solutions: Wall conditions on the sides
!           Number of iterations are for CFL of 0.5 and for
!           the rusanov solver
!           ------------------------------------------------
!
            INTEGER                            :: iterations(3:7) = [2897, 3550, 4524, 5540, 5996]
            REAL(KIND=RP), DIMENSION(3:7)      :: errors = [1.1662237969747302E-003, 3.8707028986939562E-004, &
                                                            1.0823245094648826E-004, 3.5514459858276837E-005, &
                                                            1.1953826232868892E-005]
            REAL(KIND=RP), DIMENSION(3:7)      :: residuals = [9.9114455962827790E-011, 9.9692669580629353E-011, &
                                                               9.8550101132040978E-011, 9.8967441182940477E-011, &
                                                               9.9582661331228551E-011]
!
            N = sem % spA % N
            
            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()
            
            CALL FTAssertEqual(expectedValue = iterations(N), &
                               actualValue   =  sem % numberOfTimeSteps, &
                               msg           = "Number of time steps to tolerance")
            CALL FTAssertEqual(expectedValue = residuals(N), &
                               actualValue   = sem % maxResidual, &
                               tol           = 1.d-3, &
                               msg           = "Final maximum residual")
            
            ALLOCATE(QExpected(0:sem % spA % N,0:sem % spA % N,0:sem % spA % N,N_EQN))
            
            maxError = 0.0_RP
            DO eID = 1, SIZE(sem % mesh % elements)
               DO k = 0, sem % spA % N
                  DO j = 0, sem % spA % N
                     DO i = 0, sem % spA % N 
                        CALL pointSourceFlowSolution( sem % mesh % elements(eID) % geom % x(:,i,j,k), &
                                                      QExpected(i,j,k,1:N_EQN), success )
                     END DO
                  END DO
               END DO
               maxError = MAXVAL(ABS(QExpected - sem % mesh % elements(eID) % Q))
            END DO
            CALL FTAssertEqual(expectedValue = ERRORs(N), &
                               actualValue   = maxError, &
                               tol           = 1.d-5, &
                               msg           = "Maximum error")
            
            
            CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)
   
            IF ( sharedManager % numberOfAssertionFailures() == 0 )     THEN
               WRITE(6,*) testName, " ... Passed"
            ELSE
               WRITE(6,*) testName, " ... Failed"
               WRITE(6,*) "NOTE: Failure is expected when the max eigenvalue procedure is changed."
               WRITE(6,*) "      If that is done, re-compute the expected values and modify this procedure"
            END IF 
            WRITE(6,*)
            
            CALL finalizeSharedAssertionsManager
            CALL detachSharedAssertionsManager
            
         END SUBROUTINE UserDefinedFinalize
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedTermination
!
!        -----------------------------------------------
!        Called at the the end of the main driver after 
!        everything else is done.
!        -----------------------------------------------
!
         IMPLICIT NONE  
      END SUBROUTINE UserDefinedTermination
      
      END MODULE UserDefinedFunctions
!
!//////////////////////////////////////////////////////////////////////// 
! 
   SUBROUTINE pointSourceFlowSolution(x, Q, success)
      USE UserDefinedDataStorage
      USE SMConstants
      USE Physics
      IMPLICIT NONE  
!
!           ---------
!           Arguments
!           ---------
!
      REAL(KIND=RP) :: x(3)
      REAL(KIND=RP) :: Q(*)
      LOGICAL       :: success
!
!           ---------------
!           Local variables
!           ---------------
!
      REAL(KIND=RP)                         :: tggm1, fr, rho, ff, ffp, delt
      REAL(KIND=RP)                         :: p, velocity, u, v, qq
      
      REAL(KIND=RP)                         :: tol
      CHARACTER(LEN=STRING_CONSTANT_LENGTH) :: msg
      INTEGER                               :: k
!
!     -------------------------------------------
!     Compute flow quantitites inside the element
!     The code below gives the exact solution for
!     a subsonic point source. This solution
!     requires the inflow flux, which is computed
!     in the UserDefinedFinalSetup routine at the
!     beginning of the computation.
!     -------------------------------------------
!
      success = .TRUE.
      tol   = 100.0_RP*EPSILON(1.0_RP)
      tggm1 = 2.0_RP*gamma/(gamma-1.0_RP)
      fr    = f*f/(x(1)*x(1) + x(2)*x(2))
!
!     -------------------------------------
!     Get an initial guess for the solution
!     -------------------------------------
!      
      IF(mach < 1.0_RP)     THEN
         rho = gammaMinus1/(4._RP*gamma)*(h + SQRT(h*h - 4._RP*fr))
      ELSE
         rho = 1.0_RP
         DO k = 1,10
            rho = SQRT(fr/(h - tggm1*rho**0.4_RP))
         END DO
      ENDIF
!
!     ---------------
!     Newton's method
!     ---------------
!
      DO k = 1,15
      
         ff   = tggm1*rho**gammaMinus1 + fr/rho**2 - h
         ffp  = 2._RP*gamma*rho**(gamma-2._RP) - 2._RP*fr/rho
         delt = -ff/ffp
         
         IF( abs(delt) <= tol)     EXIT
         
         rho = rho + delt
      END DO

      IF ( abs(delt) > tol )     THEN
         PRINT *, "Newton iteration on initial condition not convergedat (x,y) = ", &
                   x(1),x(2),x(3),". Delta = ",delt
         success = .FALSE.
         RETURN
      END IF
!
!     ---------------
!     Set up solution
!     ---------------
!
      Q(1) = rho
      u    = x(1)*f/rho/(x(1)*x(1) + x(2)*x(2)) 
      v    = x(2)*f/rho/(x(1)*x(1) + x(2)*x(2))
      Q(2) = rho*u
      Q(3) = rho*v
      Q(4) = 0.0_RP
      p    = rho**gamma
      Q(5) = p/gammaMinus1 + 0.5_RP*rho*(u**2 + v**2)
      
   END SUBROUTINE pointSourceFlowSolution
!
!=====================================================================================================
!=====================================================================================================
!
!
      SUBROUTINE externalStateForBoundaryName( x, t, nHat, Q, boundaryType )
!
!     ----------------------------------------------
!     Set the boundary conditions for the mesh by
!     setting the external state for each boundary.
!     ----------------------------------------------
!
      USE BoundaryConditionFunctions
      USE UserDefinedDataStorage
      USE MeshTypes
      
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: Q(N_EQN)
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryType
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP)   :: pExt
      LOGICAL         :: success
      
      IF ( boundarytype == implementedBCNames(FREE_SLIP_WALL_INDEX) )              THEN
         CALL FreeSlipWallState( x, t, nHat, Q )
      ELSE IF ( boundaryType == implementedBCNames(NO_SLIP_ADIABATIC_WALL_INDEX) ) THEN 
         CALL  NoSlipAdiabaticWallState( x, t, Q)
      ELSE IF ( boundarytype == implementedBCNames(NO_SLIP_ISOTHERMAL_WALL_INDEX)) THEN 
         CALL NoSlipIsothermalWallState( x, t, Q )
      ELSE IF ( boundaryType == implementedBCNames(OUTFLOW_SPECIFY_P_INDEX) )      THEN 
         pExt =  ExternalPressure()
         CALL ExternalPressureState ( x, t, nHat, Q, pExt )
      ELSE
         CALL pointSourceFlowSolution( x, Q, success)
      END IF

      END SUBROUTINE externalStateForBoundaryName
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ExternalGradientForBoundaryName( x, t, nHat, U_x, U_y, U_z, boundaryType )
!
!     ------------------------------------------------
!     Set the boundary conditions for the mesh by
!     setting the external gradients on each boundary.
!     ------------------------------------------------
!
      USE BoundaryConditionFunctions
      USE MeshTypes
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: U_x(N_GRAD_EQN), U_y(N_GRAD_EQN), U_z(N_GRAD_EQN)
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryType
!
!     ---------------
!     Local variables
!     ---------------
!
      IF ( boundarytype == implementedBCNames(FREE_SLIP_WALL_INDEX) )                   THEN
         CALL FreeSlipNeumann( x, t, nHat, U_x, U_y, U_z )
      ELSE IF ( boundaryType == implementedBCNames(NO_SLIP_ADIABATIC_WALL_INDEX) )       THEN 
         CALL  NoSlipAdiabaticWallNeumann( x, t, nHat, U_x, U_y, U_z)
      ELSE IF ( boundarytype == implementedBCNames(NO_SLIP_ISOTHERMAL_WALL_INDEX))       THEN 
         CALL NoSlipIsothermalWallNeumann( x, t, nHat, U_x, U_y, U_z )
      ELSE
         CALL UniformFlowNeumann( x, t, nHat, U_x, U_y, U_z )
      END IF

      END SUBROUTINE ExternalGradientForBoundaryName

