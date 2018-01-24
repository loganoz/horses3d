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
            USE HexMeshClass
            use PhysicsStorage
            IMPLICIT NONE
            CLASS(HexMesh)             :: mesh
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedInitialCondition(mesh, thermodynamics_, &
                                                      dimensionless_, &
                                                          refValues_  )
            USE SMConstants
            use PhysicsStorage
            use HexMeshClass
            implicit none
            class(HexMesh)                        :: mesh
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
!
!           ---------------
!           Local variables
!           ---------------
!
            REAL(KIND=RP) :: x(3)        
            INTEGER       :: i, j, k, eID
            REAL(KIND=RP) :: rho , u , v , w , p
            REAL(KIND=RP) :: L, u_0, rho_0, p_0
            integer       :: Nx, Ny, Nz
#if defined(NAVIERSTOKES)
            
            L     = 1.0_RP
            u_0   = 1.0_RP
            rho_0 = 1.0_RP 
            p_0   = 100.0_RP

            associate( gamma => thermodynamics_ % gamma ) 
            DO eID = 1, SIZE(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)

               DO k = 0, Nz
                  DO j = 0, Ny
                     DO i = 0, Nx 

                         x = mesh % elements(eID) % geom % x(:,i,j,k)
                       
                         rho = rho_0
                         u   =  u_0 * sin(x(1)/L) * cos(x(2)/L) * cos(x(3)/L) 
                         v   = -u_0 * cos(x(1)/L) * sin(x(2)/L) * cos(x(3)/L)
                         w   =  0.0_RP
                         p   =   p_0 + rho_0 / 16.0_RP * (                          &
                               cos(2.0_RP*x(1)/L)*cos(2.0_RP*x(3)/L) +                  &
                               2.0_RP*cos(2.0_RP*x(2)/L) + 2.0_RP*cos(2.0_RP*x(1)/L) +  &
                               cos(2.0_RP*x(2)/L)*cos(2.0_RP*x(3)/L)                    &
                               )

                         mesh % elements(eID) % storage % Q(1,i,j,k) = rho
                         mesh % elements(eID) % storage % Q(2,i,j,k) = rho*u
                         mesh % elements(eID) % storage % Q(3,i,j,k) = rho*v
                         mesh % elements(eID) % storage % Q(4,i,j,k) = rho*w
                         mesh % elements(eID) % storage % Q(5,i,j,k) = p / (gamma - 1.0_RP) + 0.5_RP * rho * (u*u + v*v + w*w)

                     END DO
                  END DO
               END DO 
               
            END DO 
            end associate
#endif            
            
         END SUBROUTINE UserDefinedInitialCondition

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
         subroutine UserDefinedSourceTerm(mesh, time, thermodynamics_, dimensionless_, refValues_)
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
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
!           Usage example
!           -------------
!           do eID = 1, mesh % no_of_elements
!              associate ( e => mesh % elements(eID) )
!              do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
!                 associate(x => e % geom % x(1,i,j,k), &
!                           y => e % geom % x(2,i,j,k), &
!                           z => e % geom % x(3,i,j,k)  )
!                 e % storage % S(:,i,j,k) = x + y + z + time
!                 end associate
!              end do                  ; end do                ; end do
!              end associate
!           end do
   
         end subroutine UserDefinedSourceTerm
!
!//////////////////////////////////////////////////////////////////////// 
! 

         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, monitors)
!
!           ----------------------------------------------------------
!           Called at the output interval to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            USE HexMeshClass
            use MonitorsClass
            IMPLICIT NONE
            CLASS(HexMesh)  :: mesh
            REAL(KIND=RP) :: time
            type(Monitor_t),  intent(in)  :: monitors
!
!           ---------------
!           Local variables
!           ---------------
!
            LOGICAL, SAVE       :: created = .FALSE.
            REAL(KIND=RP)       :: k_energy = 0.0_RP
            INTEGER             :: fUnit
            INTEGER, SAVE       :: restartfile = 0
            INTEGER             :: restartUnit
            CHARACTER(LEN=100)  :: fileName
            
         END SUBROUTINE UserDefinedPeriodicOperation
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
            use FTAssertions
            USE HexMeshClass
            use MonitorsClass
            use PhysicsStorage
            use MPI_Process_Info
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
            type(Thermodynamics_t),    intent(in) :: thermodynamics_
            type(Dimensionless_t),     intent(in) :: dimensionless_
            type(RefValues_t),         intent(in) :: refValues_
            type(Monitor_t),           intent(in) :: monitors
            real(kind=RP),          intent(in) :: elapsedTime
            real(kind=RP),          intent(in) :: CPUTime
!
!           ---------------
!           Local variables
!           ---------------
!
#if defined(NAVIERSTOKES)
            CHARACTER(LEN=29)                  :: testName           = "Taylor-Green vortex"
            REAL(KIND=RP)                      :: maxError
            REAL(KIND=RP), ALLOCATABLE         :: QExpected(:,:,:,:)
            INTEGER                            :: eID
            INTEGER                            :: i, j, k, N
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
            integer                            :: rank
            real(kind=RP), parameter           :: kinEn = 0.12499758737106952_RP
            real(kind=RP), parameter           :: kinEnRate = -4.2807169311969659E-004_RP
            real(kind=RP), parameter           :: enstrophy = 0.37499411028501956_RP 
            real(kind=RP), parameter           :: res(5) = [5.2294183691232104E-005_RP, &
                                                            0.12783424260634596_RP, &
                                                            0.12783424273963268_RP, &
                                                            0.24980299744783380_RP, &
                                                            0.61006093083852786_RP ]
            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()
            
            CALL FTAssertEqual(expectedValue = monitors % residuals % values(1,1) + 1.0_RP, &
                               actualValue   = res(1) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "continuity residual")

            CALL FTAssertEqual(expectedValue = monitors % residuals % values(2,1) + 1.0_RP, &
                               actualValue   = res(2) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "x-momentum residual")

            CALL FTAssertEqual(expectedValue = monitors % residuals % values(3,1) + 1.0_RP, &
                               actualValue   = res(3) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "y-momentum residual")

            CALL FTAssertEqual(expectedValue = monitors % residuals % values(4,1) + 1.0_RP, &
                               actualValue   = res(4) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "z-momentum residual")

            CALL FTAssertEqual(expectedValue = monitors % residuals % values(5,1) + 1.0_RP, &
                               actualValue   = res(5) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "energy residual")

            CALL FTAssertEqual(expectedValue = monitors % volumeMonitors(1) % values(1), &
                               actualValue   = kinEn, &
                               tol           = 1.0e-11_RP, &
                               msg           = "Kinetic Energy")

            CALL FTAssertEqual(expectedValue = monitors % volumeMonitors(2) % values(1) + 1.0_RP, &
                               actualValue   = kinEnRate + 1.0_RP, &
                               tol           = 1.0e-11_RP, &
                               msg           = "Kinetic Energy Rate")

            CALL FTAssertEqual(expectedValue = monitors % volumeMonitors(3) % values(1), &
                               actualValue   = enstrophy, &
                               tol           = 1.0e-11_RP, &
                               msg           = "Enstrophy")

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
      
