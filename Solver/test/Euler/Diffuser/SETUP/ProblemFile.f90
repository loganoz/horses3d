!
!////////////////////////////////////////////////////////////////////////
!
!      The Problem File contains user defined procedures
!      that are used to "personalize" i.e. define a specific
!      problem to be solved. These procedures include initial conditions,
!      exact solutions (e.g. for tests), etc. and allow modifications 
!      without having to modify the main code.
!
!      The procedures, *even if empty* that must be defined are
!
!      UserDefinedSetUp
!      UserDefinedInitialCondition(mesh)
!      UserDefinedPeriodicOperation(mesh)
!      UserDefinedFinalize(mesh)
!      UserDefinedTermination
!
!//////////////////////////////////////////////////////////////////////// 
! 
!#include "Includes.h"
      MODULE UserDefinedDataStorage
         USE SMConstants
         IMPLICIT NONE 
         REAL(KIND=RP) :: rad0, f, h 
      END MODULE UserDefinedDataStorage

         SUBROUTINE UserDefinedStartup
!
!        --------------------------------
!        Called before any other routines
!        --------------------------------
!
            use SMConstants
            IMPLICIT NONE  
!#if (!defined(NAVIERSTOKES))
!            print*, "This test case only works with NS"
!            !errorMessage(STD_OUT)
!            error stop
!#endif
         END SUBROUTINE UserDefinedStartup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalSetup(mesh &
#if defined(NAVIERSTOKES)
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#if defined(CAHNHILLIARD)
                                        , multiphase_ &
#endif
                                        )
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            USE HexMeshClass
            use UserDefinedDataStorage
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
            class(HexMesh)                      :: mesh
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           Local variables
!           ---------------
!
            INTEGER       :: nodeID
            REAL(KIND=RP) :: x(3)
            REAL(KIND=RP) :: rad

#if defined(NAVIERSTOKES)
            associate ( gammaMinus1Div2 => thermodynamics_ % gammaMinus1Div2, &
                        sqrtGamma => thermodynamics_ % sqrtGamma, &
                        gamma => thermodynamics_ % gamma, &
                        Mach => dimensionless_ % Mach )

            rad0 = HUGE(1.0_RP)
            DO nodeID = 1, SIZE(mesh % nodes)
               x   = mesh % nodes(nodeID) % x
               rad = SQRT(x(1)**2 + x(2)**2)
               rad0  = MIN(rad0, rad)
            END DO
            
            f = sqrtGamma*rad0*mach
            h = GAMMA *(1.0_RP/gammaMinus1Div2 + mach**2)

            end associate
#endif

         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         subroutine UserDefinedInitialCondition(mesh &
#if defined(NAVIERSTOKES)
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#if defined(CAHNHILLIARD)
                                        , multiphase_ &
#endif
                                        )
!
!           ------------------------------------------------
!           called to set the initial condition for the flow
!              - by default it sets an uniform initial
!                 condition.
!           ------------------------------------------------
!
            use smconstants
            use physicsstorage
            use hexmeshclass
            use fluiddata
            implicit none
            class(hexmesh)                      :: mesh
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           local variables
!           ---------------
!
            logical        :: success
            integer        :: eid, i, j, k
            real(kind=RP)  :: qq, u, v, w, p
#if defined(NAVIERSTOKES)
            real(kind=RP)  :: Q(NCONS), phi, theta
#endif

#if defined(NAVIERSTOKES) & (!defined(SPALARTALMARAS))
            DO eID = 1, SIZE(mesh % elements)
               DO k = 0, mesh % elements(eID) % Nxyz(3)
                  DO j = 0, mesh % elements(eID) % Nxyz(2)
                     DO i = 0, mesh % elements(eID) % Nxyz(1)
                        CALL pointSourceFlowSolution( mesh % elements(eID) % geom % x(:,i,j,k), &
                                                      mesh % elements(eID) % storage % Q(:,i,j,k), success, &
                                                      thermodynamics_, dimensionless_, refValues_)
                        IF(.NOT. success) error stop "Unable to compute initial condition"       
                     END DO
                  END DO
               END DO 
            END DO 
#endif

         end subroutine UserDefinedInitialCondition
#if defined(NAVIERSTOKES)
         subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
!
!           -------------------------------------------------
!           Used to define an user defined boundary condition
!           -------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: Q(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
!
!           ---------------
!           Local variables            
!           ---------------
!
            logical  :: success
#if defined(NAVIERSTOKES) & (!defined(SPALARTALMARAS))
            call pointSourceFlowSolution(x, Q, success, thermodynamics_, dimensionless_, refValues_)
#endif 
         end subroutine UserDefinedState1

         subroutine UserDefinedGradVars1(x, t, nHat, Q, U, GetGradients, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            use VariableConversion, only: GetGradientValues_f
            implicit none
            real(kind=RP), intent(in)          :: x(NDIM)
            real(kind=RP), intent(in)          :: t
            real(kind=RP), intent(in)          :: nHat(NDIM)
            real(kind=RP), intent(in)          :: Q(NCONS)
            real(kind=RP), intent(inout)       :: U(NGRAD)
            procedure(GetGradientValues_f)     :: GetGradients
            type(Thermodynamics_t), intent(in) :: thermodynamics_
            type(Dimensionless_t),  intent(in) :: dimensionless_
            type(RefValues_t),      intent(in) :: refValues_
         end subroutine UserDefinedGradVars1

         subroutine UserDefinedNeumann1(x, t, nHat, U_x, U_y, U_z)
!
!           --------------------------------------------------------
!           Used to define a Neumann user defined boundary condition
!           --------------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: U_x(NGRAD)
            real(kind=RP), intent(inout)  :: U_y(NGRAD)
            real(kind=RP), intent(inout)  :: U_z(NGRAD)
         end subroutine UserDefinedNeumann1
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, dt, Monitors)
!
!           ----------------------------------------------------------
!           Called at the output interval to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use UserDefinedDataStorage
#if defined(NAVIERSTOKES)
            use MonitorsClass
#endif
            IMPLICIT NONE
            class(HexMesh)               :: mesh
            real(kind=RP)                :: time
            real(kind=RP)                :: dt
#if defined(NAVIERSTOKES)
            type(Monitor_t), intent(in) :: monitors
#else
            logical, intent(in) :: monitors
#endif
            
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
#if defined(NAVIERSTOKES)
         subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_)
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use UserDefinedDataStorage
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
            real(kind=RP),             intent(in)  :: x(NDIM)
            real(kind=RP),             intent(in)  :: Q(NCONS)
            real(kind=RP),             intent(in)  :: time
            real(kind=RP),             intent(inout) :: S(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
!
!           Usage example
!           -------------
!           S(:) = x(1) + x(2) + x(3) + time
            S    = 0.0_RP
   
         end subroutine UserDefinedSourceTermNS
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual &
#if defined(NAVIERSTOKES)
                                                    , thermodynamics_ &
                                                    , dimensionless_  &
                                                    , refValues_ & 
#endif   
#if defined(CAHNHILLIARD)
                                                    , multiphase_ &
#endif
                                                    , monitors, &
                                                      elapsedTime, &
                                                      CPUTime   )
!
!           --------------------------------------------------------
!           Called after the solution computed to allow, for example
!           error tests to be performed
!           --------------------------------------------------------
!
            use SMConstants
            USE FTAssertions
            USE HexMeshClass
            use UserDefinedDataStorage
            use PhysicsStorage
            use FluidData
            use MonitorsClass
            IMPLICIT NONE
            class(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)    :: thermodynamics_
            type(Dimensionless_t),  intent(in)    :: dimensionless_
            type(RefValues_t),      intent(in)    :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)    :: multiphase_
#endif
            type(Monitor_t),        intent(in)    :: monitors
            real(kind=RP),             intent(in) :: elapsedTime
            real(kind=RP),             intent(in) :: CPUTime
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
#if defined(NAVIERSTOKES) &(!(SPALARTALMARAS))
            interface
               SUBROUTINE pointSourceFlowSolution(x, Q, success, thermodynamics_, &
                                                                 dimensionless_, &
                                                                     refValues_  )
                  USE UserDefinedDataStorage
                  USE SMConstants
                  USE PhysicsStorage
                  use FluidData
                  IMPLICIT NONE  
                  REAL(KIND=RP) :: x(3)
                  REAL(KIND=RP) :: Q(NCONS)
                  LOGICAL       :: success
                  type(Thermodynamics_t), intent(in)  :: thermodynamics_
                  type(Dimensionless_t),  intent(in)  :: dimensionless_
                  type(RefValues_t),      intent(in)  :: refValues_
               end subroutine pointSourceFlowSolution
            end interface
#endif
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
#if defined(NAVIERSTOKES) &(!(SPALARTALMARAS))
!
            N = mesh % elements(1) % Nxyz(1) ! This works here because all the elements have the same order
            
            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()
            
            CALL FTAssertEqual(expectedValue = iterations(N), &
                               actualValue   =  iter, &
                               msg           = "Number of time steps to tolerance")
            CALL FTAssertEqual(expectedValue = residuals(N), &
                               actualValue   = maxResidual, &
                               tol           = 1.d-3, &
                               msg           = "Final maximum residual")
            
            ALLOCATE(QExpected(NCONS,0:N,0:N,0:N))
            
            maxError = 0.0_RP
            DO eID = 1, SIZE(mesh % elements)
               DO k = 0, mesh % elements(eID) % Nxyz(3)
                  DO j = 0, mesh % elements(eID) % Nxyz(2)
                     DO i = 0, mesh % elements(eID) % Nxyz(1)
                        CALL pointSourceFlowSolution( mesh % elements(eID) % geom % x(:,i,j,k), &
                                                      QExpected(1:NCONS,i,j,k), success, &
                                                      thermodynamics_, dimensionless_, refValues_ )
                     END DO
                  END DO
               END DO
               maxError = MAXVAL(ABS(QExpected - mesh % elements(eID) % storage % Q))
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
               error stop 99
            END IF 
            WRITE(6,*)
            
            CALL finalizeSharedAssertionsManager
            CALL detachSharedAssertionsManager
#endif

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
#if defined(NAVIERSTOKES) & (!(SPALARTALMARAS))
   SUBROUTINE pointSourceFlowSolution(x, Q, success, thermodynamics_, &
                                                     dimensionless_, &
                                                         refValues_  )
      USE UserDefinedDataStorage
      USE SMConstants
      USE PhysicsStorage
      use FluidData
      IMPLICIT NONE  
!
!           ---------
!           Arguments
!           ---------
!
      REAL(KIND=RP) :: x(3)
      REAL(KIND=RP) :: Q(NCONS)
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
!     Compute flow quantities inside the element
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
#endif
!
!      