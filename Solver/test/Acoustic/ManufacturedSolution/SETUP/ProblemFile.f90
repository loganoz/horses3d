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
         subroutine UserDefinedInitialCondition(mesh)
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
!
!           ---------------
!           Local variables
!           ---------------
!
            REAL(KIND=RP) :: x(3)        
            INTEGER       :: i, j, k, eID
            real(KIND=RP) :: t ! Initial time
            REAL(KIND=RP) :: rho , u , v , w , p
            REAL(KIND=RP) :: L, u_0, rho_0, p_0
            integer       :: Nx, Ny, Nz

            integer :: myfile

            myfile = 10
            open(unit=myfile, file="output_userDefinedIC.txt", status="replace")
            
            L     = 1.0_RP
            u_0   = 1.0_RP
            rho_0 = 1.0_RP 
            p_0   = 100.0_RP

            t = 0.0_rp

            DO eID = 1, SIZE(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)

               DO k = 0, Nz
                  DO j = 0, Ny
                     DO i = 0, Nx 

                        x = mesh % elements(eID) % geom % x(:,i,j,k)

                        call evalManufacturedSolution(x(1), x(2), x(3), t, mesh % elements(eID) % storage % Q(:,i,j,k))

                        write(myfile,'(F24.16)') mesh % elements(eID) % storage % Q(:,i,j,k)

                     END DO
                  END DO
               END DO 
               
            END DO  
            
            close(myfile)

         end subroutine UserDefinedInitialCondition
#if defined(ACOUSTIC)
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
#if defined(ACOUSTIC)
         subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_, Qbase, Lambbase, Lamb_NS)
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use PhysicsStorage_CAA
            use FluidData
            IMPLICIT NONE
            real(kind=RP),             intent(in)  :: x(NDIM)
            real(kind=RP),             intent(in)  :: Q(NCONS)
            real(kind=RP),             intent(in)  :: time
            real(kind=RP),             intent(inout) :: S(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
            real(kind=RP),             intent(in)  :: Qbase(NCONSB)
            real(kind=RP),             intent(in)  :: Lambbase(NDIM)
            real(kind=RP),             intent(in)  :: Lamb_NS(NDIM)
!
!           Usage example
!           -------------
!           S(:) = x(1) + x(2) + x(3) + time

            real(rp) :: R, U, V, W, C, qm(NDIM)

            qm = Lambbase - Lamb_NS

            R = Qbase(IBRHO)
            U = Qbase(IBU)
            V = Qbase(IBV)
            W = Qbase(IBW)
            C = Qbase(IBA2)
            S(ICAARHO) = 0.0_rp
            S(ICAAU) = -qm(1) + (R*(U*sin(x(1))*sin(x(2))*sin(x(3)) + V*sin(x(3))*cos(x(1))*cos(x(2)) - W*sin(x(2))*cos(x(1))*cos(x(3)) + sin(x(2))*sin(x(3))*cos(x(1)))*exp(time) - sin(x(1)))*exp(-2.0_rp*time)/R
            S(ICAAV) = -qm(2) + (R*(-U*sin(x(3))*cos(x(1))*cos(x(2)) - V*sin(x(1))*sin(x(2))*sin(x(3)) - W*sin(x(1))*cos(x(2))*cos(x(3)) - sin(x(1))*sin(x(3))*cos(x(2)))*exp(time) - sin(x(2)))*exp(-2.0_rp*time)/R
            S(ICAAW) = -qm(3) + (R*(-U*sin(x(2))*cos(x(1))*cos(x(3)) + V*sin(x(1))*cos(x(2))*cos(x(3)) + W*sin(x(1))*sin(x(2))*sin(x(3)) + sin(x(1))*sin(x(2))*cos(x(3)))*exp(time) - sin(x(3)))*exp(-2.0_rp*time)/R
            S(ICAAP) = (C*R*exp(time)*sin(x(1))*sin(x(2))*sin(x(3)) - U*sin(x(1)) - V*sin(x(2)) - W*sin(x(3)) - 2.0_rp*cos(x(1)) - 2.0_rp*cos(x(2)) - 2.0_rp*cos(x(3)))*exp(-2.0_rp*time)
   
         end subroutine UserDefinedSourceTermNS
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual &
#ifdef ACOUSTIC
                                                    , thermodynamics_ &
                                                    , dimensionless_  &
                                                    , refValues_ &  
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
            USE NodalStorageClass               , only: NodalStorage
            IMPLICIT NONE
            class(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
#ifdef ACOUSTIC
            type(Thermodynamics_t), intent(in)    :: thermodynamics_
            type(Dimensionless_t),  intent(in)    :: dimensionless_
            type(RefValues_t),      intent(in)    :: refValues_
#endif
            type(Monitor_t),        intent(in)    :: monitors
            real(kind=RP),             intent(in) :: elapsedTime
            real(kind=RP),             intent(in) :: CPUTime
!
!           ---------------
!           Local variables
!           ---------------
!
            CHARACTER(LEN=29)                  :: testName           = "Periodic APE-4C"
            REAL(KIND=RP)                      :: maxError
            REAL(KIND=RP), ALLOCATABLE         :: QExpected(:,:,:,:)
            INTEGER                            :: eID
            INTEGER                            :: i, j, k, Nx, Ny, Nz
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
            integer                            :: rank
            ! real(kind=RP), parameter           :: kinEn = 1.2499879367819486E-01_RP
            ! real(kind=RP), parameter           :: kinEnRate = -4.2807806718622574E-04_RP
            ! real(kind=RP), parameter           :: enstrophy = 3.7499683882517909E-01_RP 
            ! real(kind=RP), parameter           :: res(5) = [1.6417830052388520E-05_RP, &
            !                                                 1.2677577061211545E-01_RP, &
            !                                                 1.2677577048633804E-01_RP, &
            !                                                 2.4981129585617484E-01_RP, &
            !                                                 6.2174425106488129E-01_RP ]
            REAL(KIND=RP) :: u(NCONS), uh(NCONS), x(NDIM), wi, wj, wk
            REAL(KIND=RP) :: error_mesh, error_elem
            
            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()

            ! AJRTODO
            print *, "HELLO"
            error_mesh = 0.0_rp
            DO eID = 1, SIZE(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)

               error_elem = 0.0_rp
               DO k = 0, Nz
                  DO j = 0, Ny
                     DO i = 0, Nx
                        x = mesh % elements(eID) % geom % x(:,i,j,k)
                        call evalManufacturedSolution(x(1),x(2),x(3),time,u)
                        uh = mesh % elements(eID) % storage % Q(:,i,j,k)
                        ! wi = NodalStorage(Nx) % w(i)
                        ! wj = NodalStorage(Ny) % w(j)
                        ! wk = NodalStorage(Nz) % w(k)
                        ! error_elem = error_elem + norm2(u - uh) * mesh % elements(eID) % geom % jacobian(i,j,k) * wi * wj * wk
                        error_elem = error_elem + norm2(u - uh)
                     END DO
                  END DO
               END DO
               error_mesh = error_mesh + error_elem / (Nx*Ny*Nz)
            END DO
            
            error_mesh = sqrt(error_mesh / SIZE(mesh % elements))
            print *, "Error in L2 norm: ", error_mesh


            ! block
            !    integer :: myfile, nxx, nyy, nzz
            !    myfile = 10
            !    print *, "FINAL TIME: ", time
            !    open(unit=myfile, file="output_MSFinal.txt", status="replace")
            !    DO eID = 1, SIZE(mesh % elements)
            !          Nxx = mesh % elements(eID) % Nxyz(1)
            !          Nyy = mesh % elements(eID) % Nxyz(2)
            !          Nzz = mesh % elements(eID) % Nxyz(3)

            !          DO k = 0, Nzz; DO j = 0, Nyy; DO i = 0, Nxx 
            !             x = mesh % elements(eID) % geom % x(:,i,j,k)
            !             call evalManufacturedSolution(x(1),x(2),x(3),time,u)
            !                               if (eID .eq. 1) then
            !          print *, "FIRST ELEMENT AT FINAL TIME (analytical)"
            !                                  print *, "coordinates"
            !          print *, x
            !          print *, "values"
            !          print *, u
            !       end if
            !             write(myfile,'(F24.16)') u
            !          end do; end do; end do;
                     
            !    END DO 
            !    close(myfile)
            !    ! error stop
            ! end block

            ! CALL FTAssertEqual(expectedValue = res(1) + 1.0_RP, &
            !                    actualValue   = monitors % residuals % values(1,1) + 1.0_RP, &
            !                    tol           = 1.0e-7_RP, &
            !                    msg           = "L2 norm")

            
            ! CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)
   
            ! IF ( sharedManager % numberOfAssertionFailures() == 0 )     THEN
            !    WRITE(6,*) testName, " ... Passed"
            !    WRITE(6,*) "This test case has no expected solution yet, only checks the residual after 5 iterations."
            ! ELSE
            !    WRITE(6,*) testName, " ... Failed"
            !    WRITE(6,*) "NOTE: Failure is expected when the max eigenvalue procedure is changed."
            !    WRITE(6,*) "      If that is done, re-compute the expected values and modify this procedure"
            !     error stop 99
            ! END IF 
            ! WRITE(6,*)
            
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

   pure subroutine evalManufacturedSolution(x, y, z, t, res)
      use SMConstants
      use PhysicsStorage_CAA
      implicit none
      real(rp), intent(in)  :: x, y, z, t
      real(rp), intent(out) :: res(NCONS)

      res(ICAARHO) = 0.0_rp
      res(ICAAU) = -cos(x)*sin(y)*sin(z)*exp(-t)
      res(ICAAV) = sin(x)*cos(y)*sin(z)*exp(-t)
      res(ICAAW) = -sin(x)*sin(y)*cos(z)*exp(-t)
      res(ICAAP) = (cos(x)+cos(y)+cos(z))*exp(-2.0_rp*t)
   end subroutine