!
!//////////////////////////////////////////////////////
!
!   @File:    ProblemFile.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Jan 23 16:27:56 2018
!   @Last revision date: Tue Feb 13 20:29:43 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 9cdffcbe5af1cc3ea1e17c83c91d73cc17fecde1
!
!//////////////////////////////////////////////////////
!
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
         SUBROUTINE UserDefinedFinalSetup(mesh , thermodynamics_, &
                                                 dimensionless_, &
                                                     refValues_ )
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            IMPLICIT NONE
            CLASS(HexMesh)                      :: mesh
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_

         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
#if defined(NAVIERSTOKES)
         subroutine userdefinedinitialcondition(mesh, thermodynamics_, &
                                                      dimensionless_, &
                                                          refvalues_  )
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
            implicit none
            class(hexmesh)                      :: mesh
            type(thermodynamics_t), intent(in)  :: thermodynamics_
            type(dimensionless_t),  intent(in)  :: dimensionless_
            type(refvalues_t),      intent(in)  :: refvalues_
!
!           ---------------
!           local variables
!           ---------------
!
            integer        :: eid, i, j, k
            real(kind=rp)  :: q(n_eqn)
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

         end subroutine userdefinedinitialcondition

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
#elif defined(CAHNHILLIARD)
         subroutine userdefinedinitialcondition(mesh, thermodynamics_, &
                                                      dimensionless_, &
                                                          refvalues_  )
            use smconstants
            use physicsstorage
            use hexmeshclass
            implicit none
            class(hexmesh)                      :: mesh
            type(thermodynamics_t), intent(in)  :: thermodynamics_
            type(dimensionless_t),  intent(in)  :: dimensionless_
            type(refvalues_t),      intent(in)  :: refvalues_
            
         end subroutine userdefinedinitialcondition
#endif
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
         end subroutine UserDefinedState1

         subroutine UserDefinedNeumann(x, t, nHat, U_x, U_y, U_z)
!
!           --------------------------------------------------------
!           Used to define a Neumann user defined boundary condition
!           --------------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: U_x(N_GRAD_EQN)
            real(kind=RP), intent(inout)  :: U_y(N_GRAD_EQN)
            real(kind=RP), intent(inout)  :: U_z(N_GRAD_EQN)
         end subroutine UserDefinedNeumann

!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, Monitors)
!
!           ----------------------------------------------------------
!           Called at the output interval to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use MonitorsClass
            IMPLICIT NONE
            CLASS(HexMesh)               :: mesh
            REAL(KIND=RP)                :: time
            type(Monitor_t), intent(in) :: monitors
            
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
         subroutine UserDefinedSourceTerm(mesh, time, thermodynamics_, dimensionless_, refValues_)
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            type(Thermodynamics_t),    intent(in) :: thermodynamics_
            type(Dimensionless_t),     intent(in) :: dimensionless_
            type(RefValues_t),         intent(in) :: refValues_
!
!           ---------------
!           Local variables
!           ---------------
!
            integer  :: i, j, k, eID
!
!           Usage example (by default no source terms are added)
!           ----------------------------------------------------
!           do eID = 1, mesh % no_of_elements
!              associate ( e => mesh % elements(eID) )
!              do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
!                 e % QDot(:,i,j,k) = e % QDot(:,i,j,k) + Source(:)
!              end do                  ; end do                ; end do
!           end do
   
         end subroutine UserDefinedSourceTerm
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual, thermodynamics_, &
                                                    dimensionless_, &
                                                        refValues_, &  
                                                          monitors, &
                                                       elapsedTime, &
                                                           CPUTime   )
!
!           --------------------------------------------------------
!           Called after the solution computed to allow, for example
!           error tests to be performed
!           --------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use FTAssertions
            use MonitorsClass
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
            type(Thermodynamics_t),    intent(in) :: thermodynamics_
            type(Dimensionless_t),     intent(in) :: dimensionless_
            type(RefValues_t),         intent(in) :: refValues_
            type(Monitor_t),          intent(in) :: monitors
            real(kind=RP),             intent(in)  :: elapsedTime
            real(kind=RP),             intent(in)  :: CPUTime
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

            CALL FTAssertEqual(expectedValue = monitors % volumeMonitors(3) % values(1), &
                               actualValue   = kinEn, &
                               tol           = 1.0e-11_RP, &
                               msg           = "Kinetic Energy")

            CALL FTAssertEqual(expectedValue = monitors % volumeMonitors(1) % values(1), &
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
                STOP 99
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
      
