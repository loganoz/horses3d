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
         SUBROUTINE UserDefinedFinalSetup(sem, thermodynamics_, &
                                                dimensionless_, &
                                                    refValues_ )
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in but before time integration
!           to allow mesh related initializations or memory allocations
!           ----------------------------------------------------------------------
!
            USE DGSEMClass
            USE Physics
            USE UserDefinedDataStorage
            IMPLICIT NONE
!
!           ---------
!           Arguments
!           ---------
!
            CLASS(DGSem)            :: sem
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
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
   
            associate ( gammaMinus1Div2 => thermodynamics_ % gammaMinus1Div2, &
                        sqrtGamma => thermodynamics_ % sqrtGamma, &
                        gamma => thermodynamics_ % gamma, &
                        Mach => dimensionless_ % Mach )

            rad0 = HUGE(1.0_RP)
            DO nodeID = 1, SIZE(sem % mesh % nodes)
               x   = sem % mesh % nodes(nodeID) % x
               rad = SQRT(x(1)**2 + x(2)**2)
               rad0  = MIN(rad0, rad)
            END DO
            
            f = sqrtGamma*rad0*mach
            h = GAMMA *(1.0_RP/gammaMinus1Div2 + mach**2)

            end associate
            
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedInitialCondition(sem, thermodynamics_, &
                                                      dimensionless_, &
                                                          refValues_ )
!
!           ------------------------------------------------
!           Called to set the initial condition for the flow
!           ------------------------------------------------
!
            USE SMConstants
            use DGSEMClass
            use PhysicsStorage
            implicit none
            class(DGSEM)                        :: sem
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
!
!           ---------------
!           Local variables
!           ---------------
!
            LOGICAL                  :: success
            INTEGER     :: i, j, k, eID
            interface
               SUBROUTINE pointSourceFlowSolution(x, Q, success, thermodynamics_, &
                                                                 dimensionless_, &
                                                                     refValues_  )
                  USE UserDefinedDataStorage
                  USE SMConstants
                  USE Physics
                  IMPLICIT NONE  
                  REAL(KIND=RP) :: x(3)
                  REAL(KIND=RP) :: Q(N_EQN)
                  LOGICAL       :: success
                  type(Thermodynamics_t), intent(in)  :: thermodynamics_
                  type(Dimensionless_t),  intent(in)  :: dimensionless_
                  type(RefValues_t),      intent(in)  :: refValues_
               end subroutine pointSourceFlowSolution
            end interface
            
            DO eID = 1, SIZE(sem % mesh % elements)
               DO k = 0, sem % mesh % elements(eID) % Nxyz(3)
                  DO j = 0, sem % mesh % elements(eID) % Nxyz(2)
                     DO i = 0, sem % mesh % elements(eID) % Nxyz(1)
                        CALL pointSourceFlowSolution( sem % mesh % elements(eID) % geom % x(:,i,j,k), &
                                                      sem % mesh % elements(eID) % storage % Q(i,j,k,1:N_EQN), success, &
                                                      thermodynamics_, dimensionless_, refValues_)
                        IF(.NOT. success) ERROR STOP "Unable to compute initial condition"       
                     END DO
                  END DO
               END DO 
            END DO 
            
         END SUBROUTINE UserDefinedInitialCondition
!
!//////////////////////////////////////////////////////////////////////// 
! 
         subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
!
!           -------------------------------------------------
!           Used to define an user defined boundary condition
!           -------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: Q(N_EQN)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
            interface
               SUBROUTINE pointSourceFlowSolution(x, Q, success, thermodynamics_, &
                                                                 dimensionless_, &
                                                                     refValues_  )
                  USE UserDefinedDataStorage
                  USE SMConstants
                  USE Physics
                  IMPLICIT NONE  
                  REAL(KIND=RP) :: x(3)
                  REAL(KIND=RP) :: Q(N_EQN)
                  LOGICAL       :: success
                  type(Thermodynamics_t), intent(in)  :: thermodynamics_
                  type(Dimensionless_t),  intent(in)  :: dimensionless_
                  type(RefValues_t),      intent(in)  :: refValues_
               end subroutine pointSourceFlowSolution
            end interface
!
!           ---------------
!           Local variables            
!           ---------------
!
            logical  :: success

            call pointSourceFlowSolution(x, Q, success, thermodynamics_, dimensionless_, refValues_)

         end subroutine UserDefinedState1

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
         SUBROUTINE UserDefinedFinalize(sem, time, thermodynamics_, dimensionless_, refValues_)
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
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
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
            interface
               SUBROUTINE pointSourceFlowSolution(x, Q, success, thermodynamics_, &
                                                                 dimensionless_, &
                                                                     refValues_  )
                  USE UserDefinedDataStorage
                  USE SMConstants
                  USE Physics
                  IMPLICIT NONE  
                  REAL(KIND=RP) :: x(3)
                  REAL(KIND=RP) :: Q(N_EQN)
                  LOGICAL       :: success
                  type(Thermodynamics_t), intent(in)  :: thermodynamics_
                  type(Dimensionless_t),  intent(in)  :: dimensionless_
                  type(RefValues_t),      intent(in)  :: refValues_
               end subroutine pointSourceFlowSolution
            end interface

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
            N = sem % mesh % elements(1) % Nxyz(1) ! This works here because all the elements have the same order
            
            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()
            
            CALL FTAssertEqual(expectedValue = iterations(N), &
                               actualValue   =  sem % numberOfTimeSteps, &
                               msg           = "Number of time steps to tolerance")
            CALL FTAssertEqual(expectedValue = residuals(N), &
                               actualValue   = sem % maxResidual, &
                               tol           = 1.d-3, &
                               msg           = "Final maximum residual")
            
            ALLOCATE(QExpected(0:N,0:N,0:N,N_EQN))
            
            maxError = 0.0_RP
            DO eID = 1, SIZE(sem % mesh % elements)
               DO k = 0, sem % mesh % elements(eID) % Nxyz(3)
                  DO j = 0, sem % mesh % elements(eID) % Nxyz(2)
                     DO i = 0, sem % mesh % elements(eID) % Nxyz(1)
                        CALL pointSourceFlowSolution( sem % mesh % elements(eID) % geom % x(:,i,j,k), &
                                                      QExpected(i,j,k,1:N_EQN), success, &
                                                      thermodynamics_, dimensionless_, refValues_ )
                     END DO
                  END DO
               END DO
               maxError = MAXVAL(ABS(QExpected - sem % mesh % elements(eID) % storage % Q))
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
               STOP 99
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
!
!//////////////////////////////////////////////////////////////////////// 
! 
   SUBROUTINE pointSourceFlowSolution(x, Q, success, thermodynamics_, &
                                                     dimensionless_, &
                                                         refValues_  )
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
      REAL(KIND=RP) :: Q(N_EQN)
      LOGICAL       :: success
      type(Thermodynamics_t), intent(in)  :: thermodynamics_
      type(Dimensionless_t),  intent(in)  :: dimensionless_
      type(RefValues_t),      intent(in)  :: refValues_
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
      associate( gammaMinus1 => thermodynamics_ % gammaMinus1, &
                 gamma => thermodynamics_ % gamma, &
                 Mach => dimensionless_ % Mach ) 

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

      end associate
      
   END SUBROUTINE pointSourceFlowSolution
!
!=====================================================================================================
!=====================================================================================================
!
!
