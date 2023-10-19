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
   module ProblemFileRoutines
      implicit none
      contains
         subroutine ReadOrderFile(filename, Nx, Ny, Nz)
            implicit none
   !
   !        ----------------------------------------------------------------------
   !        Subroutine that reads input file containing polynomial orders for mesh
   !        ----------------------------------------------------------------------
   !
            character(len=*), intent(in) :: filename          !<  Name of file containing polynomial orders to initialize
            integer, allocatable         :: Nx(:),Ny(:),Nz(:) !>  Polynomial orders for each element
            !------------------------------------------
            integer                      :: fd       ! File unit
            integer                      :: nelem    ! Number of elements
            integer                      :: i        ! counter
            !------------------------------------------
            
            open(newunit = fd, FILE = filename )   
               READ(fd,*) nelem
               
               allocate(Nx(nelem),Ny(nelem),Nz(nelem))
               
               do i = 1, nelem
                  READ(fd,*) Nx(i), Ny(i), Nz(i)
               ENDDO
            close(UNIT=fd)
            
         end subroutine ReadOrderFile
   end module
   
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
            use SMConstants
            USE HexMeshClass
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
!           Called to set the initial condition for the flow
!              - By default it sets an uniform initial
!                 condition.
!           ------------------------------------------------
!
            use smconstants
            use physicsstorage
            use hexmeshclass
            use FluidData
            implicit none
            class(HexMesh)                        :: mesh
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
            integer        :: eID, i, j, k
            real(kind=RP)  :: qq, u, v, w, p
#if defined(NAVIERSTOKES)
            real(kind=RP)  :: Q(NCONS), phi, theta
#endif

#if defined(NAVIERSTOKES)
            associate ( gammaM2 => dimensionless_ % gammaM2, &
                        gamma => thermodynamics_ % gamma )
            theta = refValues_ % AOATheta*(PI/180.0_RP)
            phi   = refValues_ % AOAPhi*(PI/180.0_RP)
      
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          Ny => mesh % elements(eID) % Nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 
                  qq = 1.0_RP
                  u  = qq*cos(theta)*COS(phi)
                  v  = qq*sin(theta)*COS(phi)
                  w  = qq*SIN(phi)
      
                  Q(1) = 1.0_RP
                  p    = 1.0_RP/(gammaM2)
                  Q(2) = Q(1)*u
                  Q(3) = Q(1)*v
                  Q(4) = Q(1)*w
                  Q(5) = p/(gamma - 1._RP) + 0.5_RP*Q(1)*(u**2 + v**2 + w**2)

                  mesh % elements(eID) % storage % Q(:,i,j,k) = Q 
               end do;        end do;        end do
               end associate
            end do

            end associate
#endif
!
!           ---------------------------------------
!           Cahn-Hilliard default initial condition
!           ---------------------------------------
!
#if defined(CAHNHILLIARD)
            call random_seed()
         
            do eid = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eid) % Nxyz(1), &
                          Ny => mesh % elements(eid) % Nxyz(2), &
                          Nz => mesh % elements(eid) % Nxyz(3) )
               associate(e => mesh % elements(eID) % storage)
               call random_number(e % c) 
               e % c = 2.0_RP * (e % c - 0.5_RP)
               end associate
               end associate
            end do
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
            use ProblemFileRoutines
            use SMConstants
            use FTAssertions
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
!
!           ---------------
!           Local variables
!           ---------------
!
            CHARACTER(LEN=29)                  :: testName           = "Euler p-adapted Cylinder"
            REAL(KIND=RP)                      :: maxError
            REAL(KIND=RP), ALLOCATABLE         :: QExpected(:,:,:,:)
            INTEGER                            :: eID
            INTEGER                            :: i, j, k
            INTEGER, allocatable               :: Nx(:), Ny(:), Nz(:)
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
!
!           -----------------------------------------------------------------------------------------
!           Expected solutions. 
!           InnerCylinder 0.0 NoSlipAdiabaticWall
!           Front 0.0 Inflow
!           bottom 0.0 FreeSlipWall
!           top 0.0 FreeSlipWall
!           Back 0.0 Inflow
!           Left 0.0 Inflow
!           Right 0.0 OutflowSpecifyP 
!           -----------------------------------------------------------------------------------------
!
!
!           ------------------------------------------------
!           Expected Solutions: Wall conditions on the sides
!           Number of iterations are for CFL of 0.3, for
!           the roe solver and mach = 0.3
!           ------------------------------------------------
!
#if defined(NAVIERSTOKES)
            real(kind=RP), parameter :: final_time = 8.9356631108400600E-002_RP
            real(kind=RP), parameter :: res(5) = [ 7.5562643059618037E-004_RP, &
                                                       7.0131047667415028E-003_RP, &
                                                       8.2577456104454747E-013_RP, &
                                                       1.0463412436767145E-002_RP, &
                                                       1.7942414864704437E-002_RP ]
            real(kind=RP), parameter :: wake_u = -9.7956966481995561E-003_RP
            real(kind=RP), parameter :: cd = 37.179086773140007_RP
            real(kind=RP), parameter :: p_aver = 7.3940194428249431_RP

            
            call ReadOrderFile('MESH/PolOrdersAfterAdaptation.omesh', Nx, Ny, Nz)
            
            success = .TRUE.
            do eID = 1, mesh % no_of_elements
               associate (e => mesh % elements(eID))
               if ( Nx(eID) /= e % Nxyz(1) .or. &
                    Ny(eID) /= e % Nxyz(2) .or. &
                    Nz(eID) /= e % Nxyz(3) ) then
                  success = .FALSE.
                  exit
               end if
               end associate
            end do
            
            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()
            
            CALL FTAssertEqual(expectedValue = final_time, &
                               actualValue   = time, &
                               tol           = 1.0e-11_RP, &
                               msg           = "Final time")

            CALL FTAssertEqual(expectedValue = res(1) + 1.0_RP, &
                               actualValue   = monitors % residuals % values(1,1) + 1.0_RP, &
                               tol           = 1.0e-11_RP, &
                               msg           = "continuity residual")

            CALL FTAssertEqual(expectedValue = res(2) + 1.0_RP, &
                               actualValue   = monitors % residuals % values(2,1) + 1.0_RP, &
                               tol           = 1.0e-11_RP, &
                               msg           = "x-momentum residual")

            CALL FTAssertEqual(expectedValue = res(3) + 1.0_RP, &
                               actualValue   = monitors % residuals % values(3,1) + 1.0_RP, &
                               tol           = 1.0e-11_RP, &
                               msg           = "y-momentum residual")

            CALL FTAssertEqual(expectedValue = res(4) + 1.0_RP, &
                               actualValue   = monitors % residuals % values(4,1) + 1.0_RP, &
                               tol           = 1.0e-11_RP, &
                               msg           = "z-momentum residual")

            CALL FTAssertEqual(expectedValue = res(5) + 1.0_RP, &
                               actualValue   = monitors % residuals % values(5,1) + 1.0_RP, &
                               tol           = 1.0e-11_RP, &
                               msg           = "energy residual")

            CALL FTAssertEqual(expectedValue = wake_u + 1.0_RP, &
                               actualValue   = monitors % probes(1) % values(1) + 1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Wake final x-velocity at the point [0,2.0,4.0]")

            CALL FTAssertEqual(expectedValue = cd, &
                               actualValue   = monitors % surfaceMonitors(1) % values(1), &
                               tol           = 1.d-10, &
                               msg           = "Drag coefficient")

            CALL FTAssertEqual(expectedValue = p_aver, &
                               actualValue   = monitors % surfaceMonitors(2) % values(1), &
                               tol           = 1.d-11, &
                               msg           = "Pressure average over the cylinder")
            
            CALL FTAssertEqual(expectedValue = .TRUE. , &
                               actualValue   = success, &
                               msg           = "Polynomial orders 'MESH/PolOrdersAfterAdaptation.omesh'")


            CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)
   
            IF ( sharedManager % numberOfAssertionFailures() == 0 )     THEN
               WRITE(6,*) testName, " ... Passed"
               WRITE(6,*) "This test case has no expected solution yet, only checks the residual after 27 iterations."
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
      