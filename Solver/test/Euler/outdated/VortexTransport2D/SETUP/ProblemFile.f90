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
#if defined(NAVIERSTOKES)
            integer        :: eid, i, j, k
            real(kind=rp)  :: q(NCONS)
            real(kind=RP)            :: r2 , rho , u , w , T, p
            real(kind=RP)            :: x(NDIM)
            interface
               subroutine VortexTransportInitialCondition(x,gamma, Mach, cv, q)
                  use SMConstants
                  implicit none
                  real(kind=RP), intent(in)  :: x(3), gamma, Mach, cv
                  real(kind=RP), intent(out) :: q(5)
               end subroutine VortexTransportInitialCondition
            end interface

            associate ( gamma => Thermodynamics_ % Gamma , Mach => Dimensionless_ % Mach , cv => Dimensionless_ % cv )

            do eid = 1, mesh % no_of_elements
               associate( nx => mesh % elements(eid) % nxyz(1), &
                          ny => mesh % elements(eid) % nxyz(2), &
                          nz => mesh % elements(eid) % nxyz(3) )
               do k = 0, nz;  do j = 0, ny;  do i = 0, nx 
                  x = mesh % elements(eID) % geom % x(:,i,j,k)
                  call VortexTransportInitialCondition(x, gamma, Mach, cv, q)
                  mesh % elements(eid) % storage % q(:,i,j,k) = q 
               end do;        end do;        end do
               end associate
            end do
            end associate
#endif
         end subroutine UserDefinedInitialCondition
#if defined(NAVIERSTOKES)
         subroutine VortexTransportInitialCondition(x,gamma, Mach, cv, q)
            use SMConstants
            implicit none
            real(kind=RP), intent(in)  :: x(3), gamma, Mach, cv
            real(kind=RP), intent(out) :: q(5)
!
!           ---------------
!           Local variables
!           ---------------
!
            real(kind=RP), parameter :: XC = 0.0_RP, ZC = 0.0_RP
            real(kind=RP), parameter :: beta = 0.05_RP, R = 0.1_RP
            real(kind=RP), parameter :: AngleOfAttack = 0.0_RP
            real(kind=RP)            :: r2, rho, u, w, T, p
            
            r2 = ((x(1) - XC)*(x(1) - XC) + (x(3) - ZC)*(x(3) - ZC)) / (R*R)
         
            u =  (cos(AngleOfAttack) - Beta * (x(3) - ZC) / R * exp(-0.5_RP * r2))
            w =  (sin(AngleOfAttack) + Beta * (x(1) - XC) / R * exp(-0.5_RP * r2))
             
            T = (1.0_RP - gamma * Mach * Mach * beta * beta / (2.0_RP * gamma * cv) * exp(-r2) )
            rho = T**cv
            p  = rho * T / (gamma * Mach * Mach)
         
            q(1) = rho
            q(2) = rho * u
            q(3) = 0.0_RP
            q(4) = rho * w
            q(5) = p/(gamma - 1.0_RP) + 0.5_RP * q(1)*(u**2 + w**2)

         end subroutine VortexTransportInitialCondition
#endif

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
            use FTAssertions
            use SMConstants
            USE HexMeshClass
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
#if defined(NAVIERSTOKES)
!
!           ---------------
!           Local variables
!           ---------------
!
            CHARACTER(LEN=25) :: testName = "Vortex transport 2D KEPEC"
            integer        :: eID, i, j, k, eq
            real(kind=RP)  :: error(5), localError(5), x(3), q(5)
            real(kind=RP), parameter   :: errorS(5) = [1.6621048545506500E-004_RP, &
                                                       4.3570834365793676E-003_RP, &
                                                       1.1681607090392261E-013_RP, &
                                                       4.5003665864759282E-003_RP, &
                                                       4.6163433923638308E-003_RP    ]
            integer, parameter         :: finalIter = 3724
            real(kind=RP), parameter   :: residuals(5) = [5.0761491133344272E-003_RP, &
                                                          0.17016555967142113_RP, &
                                                          6.7716306721775754E-014_RP, &
                                                          0.47051042285525568_RP, &
                                                          0.18049802097617726_RP ]
            real(kind=RP), parameter   :: finalTime = 10.0_RP
            real(kind=RP), parameter   :: entropy = 0.68517396155684973_RP
            real(kind=RP), parameter   :: kinEn = 0.50000805172162699_RP
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            interface
               subroutine VortexTransportInitialCondition(x,gamma, Mach, cv, q)
                  use SMConstants
                  implicit none
                  real(kind=RP), intent(in)  :: x(3), gamma, Mach, cv
                  real(kind=RP), intent(out) :: q(5)
               end subroutine VortexTransportInitialCondition
            end interface

            associate ( gamma => Thermodynamics_ % Gamma , Mach => Dimensionless_ % Mach , cv => Dimensionless_ % cv )

            error = 0.0_RP
            do eID = 1, mesh % no_of_elements
               associate(e => mesh % elements(eID))
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  x = mesh % elements(eID) % geom % x(:,i,j,k)
                  call VortexTransportInitialCondition(x, gamma, Mach, cv, q)
      
                  localError = abs(q - e % storage % Q(:,i,j,k))

                  do eq = 1, 5
                     error(eq) = max(localError(eq), error(eq))
                  end do
               end do                  ; end do                ; end do
               end associate
            end do
            end associate
   
            call initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()

            CALL FTAssertEqual(expectedValue = finalIter, &
                               actualValue   = iter, &
                               msg           = "final iter")

            CALL FTAssertEqual(expectedValue = time, &
                               actualValue   = finalTime, &
                               tol           = 1.0e-12_RP, &
                               msg           = "final time")

            CALL FTAssertEqual(expectedValue = monitors % residuals % values(1,1), &
                               actualValue   = residuals(1), &
                               tol           = 1.0e-7_RP, &
                               msg           = "continuity residual")

            CALL FTAssertEqual(expectedValue = monitors % residuals % values(2,1), &
                               actualValue   = residuals(2), &
                               tol           = 1.0e-7_RP, &
                               msg           = "x-momentum residual")

            CALL FTAssertEqual(expectedValue = monitors % residuals % values(3,1), &
                               actualValue   = residuals(3), &
                               tol           = 1.0e-7_RP, &
                               msg           = "y-momentum residual")

            CALL FTAssertEqual(expectedValue = monitors % residuals % values(4,1), &
                               actualValue   = residuals(4), &
                               tol           = 1.0e-7_RP, &
                               msg           = "z-momentum residual")

            CALL FTAssertEqual(expectedValue = monitors % residuals % values(5,1), &
                               actualValue   = residuals(5), &
                               tol           = 1.0e-7_RP, &
                               msg           = "energy residual")

            CALL FTAssertEqual(expectedValue = error(1), &
                               actualValue   = errorS(1), &
                               tol           = 1.0e-7_RP, &
                               msg           = "continuity error")

            CALL FTAssertEqual(expectedValue = error(2), &
                               actualValue   = errorS(2), &
                               tol           = 1.0e-7_RP, &
                               msg           = "x-momentum error")

            CALL FTAssertEqual(expectedValue = error(3), &
                               actualValue   = errorS(3), &
                               tol           = 1.0e-7_RP, &
                               msg           = "y-momentum error")

            CALL FTAssertEqual(expectedValue = error(4), &
                               actualValue   = errorS(4), &
                               tol           = 1.0e-7_RP, &
                               msg           = "z-momentum error")

            CALL FTAssertEqual(expectedValue = error(5), &
                               actualValue   = errorS(5), &
                               tol           = 1.0e-7_RP, &
                               msg           = "energy error")

            CALL FTAssertEqual(expectedValue = monitors % volumeMonitors(3) % values(1,1), &
                               actualValue   = kinEn, &
                               tol           = 1.0e-11_RP, &
                               msg           = "Kinetic Energy")

            CALL FTAssertEqual(expectedValue = monitors % volumeMonitors(1) % values(1,1), &
                               actualValue   = entropy, &
                               tol           = 1.0e-11_RP, &
                               msg           = "Entropy")

            CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)
   
            IF ( sharedManager % numberOfAssertionFailures() == 0 )     THEN
               WRITE(6,*) testName, " ... Passed"
               WRITE(6,*) "This test case has no expected solution yet, only checks the residual after 100 iterations."
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
      