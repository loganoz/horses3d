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
            integer        :: eid, i, j, k
            real(kind=RP)  :: qq, u, v, w, p
#if defined(NAVIERSTOKES)
            real(kind=RP)  :: Q(NCONS), phi, theta
#endif
            real(kind=RP)  :: x(3) 

!
!           ---------------------------------------
!           Navier-Stokes default initial condition
!           ---------------------------------------
!
#if defined(NAVIERSTOKES)
            associate ( gammaM2 => dimensionless_ % gammaM2, &
                        gamma => thermodynamics_ % gamma )
            theta = refvalues_ % AOAtheta*(pi/180.0_RP)
            phi   = refvalues_ % AOAphi*(pi/180.0_RP)
      
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 

                  x = mesh % elements(eID) % geom % x(:, i, j, k)

                  call random_number(u)
                  call random_number(v)
                  call random_number(w)

                  v = 0.0_RP !* (1.0_RP - x(3)**2) !+ ( v - 0.5_RP) * (1.0_RP - x(3)**2) !* (1 + (v*2 - 1) * 0.6_RP)
                  !u = 1.0_RP
                  u = 0.0_RP !( u - 0.5_RP )    !v * (u*2 - 1)*0.3_RP
                  w = 0.0_RP !( w - 0.5_RP )    !v * (w*2 - 1)*0.3_RP
      
                  q(1) = 1.0_RP
                  p    = 1.0_RP/(gammaM2)
                  q(2) = q(1)*u
                  q(3) = q(1)*v
                  q(4) = q(1)*w
                  q(5) = p/(gamma - 1._RP) + 0.5_RP*q(1)*(u**2 + v**2 + w**2)

                  mesh % elements(eID) % storage % q(:,i,j,k) = q 
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

!
!           ---------------
!           local variables
!           ---------------
!
            integer        :: eid, i, j, k
            real(kind=RP)  :: qq, u, v, w, p
#if defined(NAVIERSTOKES)
            real(kind=RP)  :: phi, theta
#endif

!
!           ---------------------------------------
!           Navier-Stokes default initial condition
!           ---------------------------------------
!
#if defined(NAVIERSTOKES)
            associate ( gammaM2 => dimensionless_ % gammaM2, &
                        gamma => thermodynamics_ % gamma )
            theta = refvalues_ % AOAtheta*(pi/180.0_RP)
            phi   = refvalues_ % AOAphi*(pi/180.0_RP)
      
            ! do eID = 1, mesh % no_of_elements
            !    associate( Nx => mesh % elements(eID) % Nxyz(1), &
            !               ny => mesh % elemeNts(eID) % nxyz(2), &
            !               Nz => mesh % elements(eID) % Nxyz(3) )
            !    do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 

            !       x = mesh % elements(eID) % geom % x(:, i, j, k)

                  call random_number(u)
                  call random_number(v)
                  call random_number(w)

                  v = 0.0_RP !* (1.0_RP - x(3)**2) !+ ( v - 0.5_RP) * (1.0_RP - x(3)**2) !* (1 + (v*2 - 1) * 0.6_RP)
                  !u = 1.0_RP
                  u = 0.0_RP !( u - 0.5_RP )    !v * (u*2 - 1)*0.3_RP
                  w = 0.0_RP !( w - 0.5_RP )    !v * (w*2 - 1)*0.3_RP
      
                  q(1) = 1.0_RP
                  p    = 1.0_RP/(gammaM2)
                  q(2) = q(1)*u
                  q(3) = q(1)*v
                  q(4) = q(1)*w
                  q(5) = p/(gamma - 1._RP) + 0.5_RP*q(1)*(u**2 + v**2 + w**2)

            !      mesh % elements(eID) % storage % q(:,i,j,k) = q 
            !    end do;        end do;        end do
            !    end associate
            ! end do

            end associate
#endif            
         end subroutine UserDefinedState1

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
!           ---------------
!           Local variables
!           ---------------
!
            integer  :: i, j, k, eID
!
!           Usage example
!           -------------
!           S(:) = x(1) + x(2) + x(3) + time
   
            S(1) = 0.0_RP
            S(2) = 0.0_RP
            S(3) = 0.0_RP !0.000355556_RP 
            S(4) = 0.0_RP
            S(5) = 0.0_RP !0.000355556_RP * Q(3) / Q(1)

            ! Reynolds centerline U=1 1012.5  -> 0.001975_RP
            ! Reynolds centerline U=1 3789.47 -> 0.000527778_RP
            ! Reynolds centerline U=1 5625.0  -> 0.000355556_RP
            ! Critical Reynolds is 5772. At this value, any perturbation should finish in turbulent flow. 
            ! at lower values, turbulence should be maintained down to lower values if transition is obtained.
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
            CHARACTER(LEN=38)                  :: testName           = "Rebound and periodicity for particles."
            REAL(KIND=RP)                      :: maxError
            REAL(KIND=RP), ALLOCATABLE         :: QExpected(:,:,:,:)
            INTEGER                            :: eID
            INTEGER                            :: i, j, k, N
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
#if defined(NAVIERSTOKES)
            REAL(KIND=RP)                      :: residuals       = 621.847759506997_RP !It uses random functions so I guess it depends on the compiler. 
            !The coded value if for Alderaan Intel release. Bender intel 2015 release gives 697.169619289106_RP

            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()

            CALL FTAssertEqual(expectedValue = residuals, &
                               actualValue   = maxResidual, &
                               tol           = 1.d-11, &
                               msg           = "Final maximum residual")


            CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)
   
            IF ( sharedManager % numberOfAssertionFailures() == 0 )     THEN
               WRITE(6,*) testName, " ... Passed"
               WRITE(6,*) "This test case checks the residual after 100 iterations."
            ELSE
               WRITE(6,*) testName, " ... Failed"
               WRITE(6,*) "NOTE: Failure is expected if particle model is modified."
               WRITE(6,*) "      If that is done, re-compute the expected values and modify this procedure."   
               WRITE(6,*) "NOTE: Failure is expected if compiler version or architecture changes."
               WRITE(6,*) "      This test case uses random functions for the injection of the particles."
               WRITE(6,*) "      The coded residual is for Alderaan Intel Release configuration."                 
               WRITE(6,*) "If it fails, the particles are supposed to be here:"
               WRITE(6,*)
               WRITE(6,*) "i,x coord,y coord,z coord,u,v,w,T,Row ID"
               WRITE(6,*) " 1,0.378817E-01,0.122848E+00,0.333146E-01,0.187350E+01,0.936817E+01,-.187351E+01,0.187460E+01"
               WRITE(6,*) " 2,0.273126E-01,0.121599E+00,0.228846E-01,0.187486E+01,0.937480E+01,-.187483E+01,0.187592E+01"
               WRITE(6,*) " 3,0.185284E-01,0.120365E+00,0.357924E-01,0.187619E+01,0.938124E+01,-.187608E+01,0.187724E+01"
               WRITE(6,*) " 4,0.222413E-01,0.119143E+00,0.177665E-01,0.187739E+01,0.938757E+01,-.187730E+01,0.187858E+01"
               WRITE(6,*) " 5,0.362321E-01,0.117931E+00,0.309296E-01,0.187861E+01,0.939416E+01,-.187875E+01,0.188013E+01"
               WRITE(6,*) " 6,0.228732E-01,0.116718E+00,0.236805E-01,0.187962E+01,0.940013E+01,-.187954E+01,0.188137E+01"
               WRITE(6,*) " 7,0.307956E-01,0.115512E+00,0.295983E-01,0.188101E+01,0.940604E+01,-.188106E+01,0.188219E+01"
               WRITE(6,*) " 8,0.187440E-01,0.114308E+00,0.274227E-01,0.188209E+01,0.941237E+01,-.188185E+01,0.188362E+01"
               WRITE(6,*) " 9,0.272797E-01,0.113103E+00,0.225712E-01,0.188405E+01,0.941854E+01,-.188370E+01,0.188499E+01"
               WRITE(6,*) "10,0.211706E-01,0.113104E+00,0.375182E-01,0.188373E+01,0.941872E+01,0.188305E+01,0.188478E+01"
               WRITE(6,*)
               WRITE(6,*) "Compare with RESULTS/Pouransari0041/Pouransari0041.parts.0000000100.csv"
               WRITE(6,*) "NOTE: If you run in other machine or with  other compiler these results may change."
               WRITE(6,*) "      In that case, visualize in paraview if particles are interacting with BC the way"
               WRITE(6,*) "      they should."
               WRITE(6,*) " Bender intel 2015 residual 697.169619289106_RP"
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
      
