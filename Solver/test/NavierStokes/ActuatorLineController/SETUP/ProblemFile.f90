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
!
!//////////////////////////////////////////////////////////////////////// 
!
Module oscarAbl  !
    use smconstants
    Implicit None
!
    private
    public initialize_abl, get_velocity_at_height, u_star, LinearInterpolate, uStar, u_zRef, kappa, z0, d, zRef, timeArray, windSpeedArray
!

    real(kind=RP) :: z ! normal distance to wall
    real(kind=RP) :: uStar ! average friction velocity
    real(kind=RP) :: kappa ! von karman constant
    real(kind=RP) :: z0 ! roughness friction velocity
    real(kind=RP) :: zRef ! reference height for ABL
    real(kind=RP) :: u_zRef ! wind speed at reference height
    real(kind=RP) :: d ! plane displacement
    real(kind=RP) :: z_change, z_final, eps, dudz0, m, b ! for zero gradient blend
    REAL(kind=RP), DIMENSION(6) :: timeArray = (/0.0, 5.0, 30.0, 50.0, 70.0, 90.0/) ! Example time values
    REAL(kind=RP), DIMENSION(6) :: windSpeedArray = (/3.0, 3.0, 15.0, 15.0, 13.0, 13.0/) ! Example wind speed values
    contains
!
    Subroutine initialize_abl()
        Implicit None
    
        ! this must match the control file definitons
        kappa    = 0.41_RP
        zRef     = 119.0_RP
        z0       = 0.0002_RP ! wall roughness for sea or lakes according to "BCA Green Mark for Residential Buildings: Technical Guide and Requirements. 2016." Available at https://www1.bca.gov.sg/docs/default-source/docs-corp-buildsg/sustainability/gm-rb-2016-technical-guide-_rev010120.pdf
        d        = 0.0_RP       ! wall displacement
        
        u_zRef = LinearInterpolate(0.0_RP, timeArray, windSpeedArray) ! wind speed at reference height for time = 0.0
        uStar = u_star(u_zRef, kappa, zRef, z0, d)        

        z_change = 325.00_RP
        eps = 0.001_RP
        z_final = 386.45_RP ! zmax domain
        dudz0 = (u_abl(z_change, uStar, kappa, d, z0) - u_abl(z_change-eps, uStar, kappa, d, z0)) / (eps)
        m = (0.0_RP-dudz0) / (z_final-z_change)
        b = -z_final*m
    End Subroutine initialize_abl
!
    Function u_star(u_zRef, kappa, zRef, z0, d) result(uStar)
        Implicit None
        real(kind=RP) :: u_zRef ! wind speed at reference height
        real(kind=RP) :: kappa ! von karman constant
        real(kind=RP) :: z0 ! roughness friction velocity
        real(kind=RP) :: zRef ! reference height for ABL
        real(kind=RP) :: uStar
        real(kind=RP) :: d ! plane displacement

        uStar = (u_zRef*kappa) / log( (zRef - d + z0) / z0 )        
    End Function u_star
!    
    Function u_abl(z,uStar,kappa,d,z0) result(u_log)
        Implicit None
        real(kind=RP) :: z ! normal distance to wall
        real(kind=RP) :: uStar ! average friction velocity
        real(kind=RP) :: kappa ! von karman constant
        real(kind=RP) :: z0 ! roughness friction velocity
        real(kind=RP) :: d ! plane displacement
        real(kind=RP) :: u_log

        if (z .gt. d) then
            u_log = uStar / kappa * log( (z - d + z0) / z0 )
        else
            u_log = 0.0_RP
        end if

    End Function u_abl
!
    Function get_velocity_at_height(z) result(u_z)
        Implicit None
            real(kind=RP), intent(in) :: z
            real(kind=RP)             :: u_z

            ! log law for the whole domain expect the upper part
            if (z .le. z_change) then
                u_z = u_abl(z, uStar, kappa, d, z0)
            else
                ! blend smooth with non zero gradient
                u_z = u_abl(z_change, uStar, kappa, d, z0) + 0.5_RP*m*(z**2-z_change**2) + b*(z-z_change)
            end if 
    End Function get_velocity_at_height

    FUNCTION LinearInterpolate(x, xArray, yArray) RESULT(y)
        REAL(kind=RP), INTENT(IN) :: x                ! The value to interpolate
        REAL(kind=RP), DIMENSION(:), INTENT(IN) :: xArray  ! Array of x values (e.g., time values)
        REAL(kind=RP), DIMENSION(:), INTENT(IN) :: yArray  ! Array of y values (e.g., wind speed values)
        REAL(kind=RP) :: y                            ! Interpolated y value
        INTEGER :: i
    
        ! Ensure x is within the range of xArray
        IF (x <= xArray(1)) THEN
           y = yArray(1)
           RETURN
        ELSE IF (x >= xArray(SIZE(xArray))) THEN
           y = yArray(SIZE(yArray))
           RETURN
        END IF
    
        ! Find the interval [xArray(i), xArray(i+1)] where x lies
        DO i = 1, SIZE(xArray) - 1
           IF (x >= xArray(i) .AND. x <= xArray(i+1)) THEN
              y = yArray(i) + (x - xArray(i)) * (yArray(i+1) - yArray(i)) / (xArray(i+1) - xArray(i))
              RETURN
           END IF
        END DO
    END FUNCTION LinearInterpolate
!
End Module oscarAbl
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
            use oscarAbl
            IMPLICIT NONE
            CLASS(HexMesh)                      :: mesh
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
#if defined(NAVIERSTOKES)
!
!           Initialize Channel variables
!           ----------------------------
            call initialize_abl()
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
            use oscarAbl
#if defined(NAVIERSTOKES)
#endif
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
!           Local variables
!           ---------------
!
#if defined(NAVIERSTOKES)
            real(kind=RP) :: x(3)        
            INTEGER       :: i, j, k, eID
            real(kind=RP) :: rho , u , v , w , p
            real(kind=RP) :: rho_0, p_0, uDim, zDim
            real(kind=RP) :: Lr ! reference length
            integer       :: Nx, Ny, Nz
            
            rho_0 = 1.0_RP 
            p_0   = 1.0_RP
            Lr    = 1.0_RP  ! in m, must match with the actual used in the control file

            associate( gamma   => thermodynamics_ % gamma , &
                       gammaM2 => dimensionless_ % gammaM2, &
                       uref    => refValues_ % V) 
            do eID = 1, size(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)

               do k = 0, Nz   ;  do j = 0, Ny   ;  do i = 0, Nx 
                  
                  x = mesh % elements(eID) % geom % x(:,i,j,k)

                  ! Warining: I assume here that the normal direction is the z one, ie x(3)
                  zDim = x(3) * Lr
                  uDim = get_velocity_at_height(zDim)

                  ! neutral abl, constant temperature, pressure and density
                  rho = rho_0
                  u   =  uDim / uref
                  v   =  0.0_RP
                  w   =  0.0_RP
                  p   =  p_0 / ( gammaM2 ) 

                  mesh % elements(eID) % storage % Q(1,i,j,k) = rho
                  mesh % elements(eID) % storage % Q(2,i,j,k) = rho*u
                  mesh % elements(eID) % storage % Q(3,i,j,k) = rho*v
                  mesh % elements(eID) % storage % Q(4,i,j,k) = rho*w
                  mesh % elements(eID) % storage % Q(5,i,j,k) = p / (gamma - 1.0_RP) + 0.5_RP * rho * (u*u + v*v + w*w)
                  
               end do         ;  end do         ;  end do
               
            end do
            end associate

            contains

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
            use oscarAbl
#if defined(_HAS_MPI_)
            use mpi
#endif

            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: Q(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_

            real(kind=RP) :: rho , u , v , w , p
            real(kind=RP) :: rho_0, p_0, uDim, zDim, tDim
            real(kind=RP) :: Lr ! reference length
#if defined(_HAS_MPI_)
            integer        :: irank, ierr
#endif
            
            associate( gamma   => thermodynamics_ % gamma , &
                       gammaM2 => dimensionless_ % gammaM2, &
                       uref    => refValues_ % V) 
                rho_0 = 1.0_RP 
                p_0   = 1.0_RP
                Lr    = 1.0_RP  ! in m, must match with the actual used in the control file
                ! Warining: I assume here that the normal direction is the z one, ie x(3)
                zDim = x(3) * Lr
                tDim = t * Lr / uref
                u_zRef = LinearInterpolate(tDim, timeArray, windSpeedArray) ! wind speed at reference height for time = t
                uStar = u_star(u_zRef, kappa, zRef, z0, d)
                uDim = get_velocity_at_height(zDim)

                ! neutral abl, constant temperature, pressure and density
                rho = rho_0
                u   =  uDim / uref
                v   =  0.0_RP
                w   =  0.0_RP
                p   =  p_0 / ( gammaM2 ) 

                Q(1) = rho
                Q(2) = rho*u
                Q(3) = rho*v
                Q(4) = rho*w
                Q(5) = p / (gamma - 1.0_RP) + 0.5_RP * rho * (u*u + v*v + w*w)
           end associate

! #if defined(_HAS_MPI_)
!                 call mpi_comm_rank(MPI_COMM_WORLD, irank, ierr)
!                 if (irank == 0) then
!                     write(STD_OUT, '(A,F7.3)') "Wind speed at z_ref: ", u_zRef
!                 end if
! #endif

         end subroutine UserDefinedState1

         subroutine UserDefinedGradVars1(x, t, nHat, Q, U, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)          :: x(NDIM)
            real(kind=RP), intent(in)          :: t
            real(kind=RP), intent(in)          :: nHat(NDIM)
            real(kind=RP), intent(in)          :: Q(NCONS)
            real(kind=RP), intent(inout)       :: U(NGRAD)
            type(Thermodynamics_t), intent(in) :: thermodynamics_
            type(Dimensionless_t),  intent(in) :: dimensionless_
            type(RefValues_t),      intent(in) :: refValues_
         end subroutine UserDefinedGradVars1

         subroutine UserDefinedNeumann1(x, t, nHat, Q, U_x, U_y, U_z, flux, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)    :: x(NDIM)
            real(kind=RP), intent(in)    :: t
            real(kind=RP), intent(in)    :: nHat(NDIM)
            real(kind=RP), intent(in)    :: Q(NCONS)
            real(kind=RP), intent(in)    :: U_x(NGRAD)
            real(kind=RP), intent(in)    :: U_y(NGRAD)
            real(kind=RP), intent(in)    :: U_z(NGRAD)
            real(kind=RP), intent(inout) :: flux(NCONS)
            type(Thermodynamics_t), intent(in) :: thermodynamics_
            type(Dimensionless_t),  intent(in) :: dimensionless_
            type(RefValues_t),      intent(in) :: refValues_
         end subroutine UserDefinedNeumann1
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, dt, Monitors)
!
!           -------------------------------------------------------------------------
!           Called before every time-step to allow periodic operations to be performed
!           Here:
!              * Reading mean momentum/velocity in the volume
!              * Obtaining fTurbulentChannel according to it
!           -------------------------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
#if defined(NAVIERSTOKES)
            use MonitorsClass
#endif
            IMPLICIT NONE
            !-arguments---------------------------------------------------
            CLASS(HexMesh)               :: mesh
            real(kind=RP)                :: time
            real(kind=RP)                :: dt
#if defined(NAVIERSTOKES)
            type(Monitor_t), intent(in) :: monitors
#else
            logical, intent(in) :: monitors
#endif
#if defined(NAVIERSTOKES)
#endif
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
#if defined(NAVIERSTOKES)
         subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_) ! , Q, dt
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
            !-arguments--------------------------------------------------
            real(kind=RP),             intent(in)  :: x(NDIM)
            real(kind=RP),             intent(in)  :: Q(NCONS)
            real(kind=RP),             intent(in)  :: time
            real(kind=RP),             intent(out) :: S(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
            
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
            CLASS(HexMesh)                        :: mesh
            real(kind=RP)                         :: time
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
            CHARACTER(LEN=29)                  :: testName           = "Actuator Line Controller"
            REAL(KIND=RP)                      :: maxError
            REAL(KIND=RP), ALLOCATABLE         :: QExpected(:,:,:,:)
            INTEGER                            :: eID
            INTEGER                            :: i, j, k, N
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
!
!           --------------------------------------------------
!           Expected Solutions: Volume integral of source term
!           --------------------------------------------------
!
#if defined(NAVIERSTOKES)
  
            real(kind=RP), parameter :: residuals(5) = [ 1.0080719192306057E-04_RP,&
                                                         6.4849063073696856E-03_RP,&
                                                         1.3641651819972784E-03_RP,&
                                                         3.1145069351486703E-03_RP,&
                                                         2.2527522885979795E-01_RP]
        

            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()

            CALL FTAssertEqual(expectedValue = residuals(1)+1.0_RP, &
                               actualValue   = monitors % residuals % values(1,1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Continuity residual")

            CALL FTAssertEqual(expectedValue = residuals(2)+1.0_RP, &
                               actualValue   = monitors % residuals % values(2,1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "X-Momentum residual")

            CALL FTAssertEqual(expectedValue = residuals(3)+1.0_RP, &
                               actualValue   = monitors % residuals % values(3,1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Y-Momentum residual")

            CALL FTAssertEqual(expectedValue = residuals(4)+1.0_RP, &
                               actualValue   = monitors % residuals % values(4,1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Z-Momentum residual")

            CALL FTAssertEqual(expectedValue = residuals(5)+1.0_RP, &
                               actualValue   = monitors % residuals % values(5,1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Energy residual")

            CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)
   
            IF ( sharedManager % numberOfAssertionFailures() == 0 )     THEN
               WRITE(6,*) testName, " ... Passed"
               WRITE(6,*) "This test case has no expected solution yet, only checks the residual after 100 iterations."
            ELSE
               WRITE(6,*) testName, " ... Failed"
               WRITE(6,*) "NOTE: This test is dependant on the external controller used."
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
      
