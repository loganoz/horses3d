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
#include "Includes.h"
#if defined(NAVIERSTOKES)
#define NNS NCONS
#define NGRADNS NGRAD
#elif defined(INCNS)
#define NNS NCONS
#define NGRADNS NCONS
#endif
module ProblemFileFunctions
   implicit none

   abstract interface
      subroutine UserDefinedStartup_f
      end subroutine UserDefinedStartup_f
   
      SUBROUTINE UserDefinedFinalSetup_f(mesh &
#if defined(NAVIERSTOKES) || defined(INCNS)
                                     , thermodynamics_ &
                                     , dimensionless_  &
                                     , refValues_ & 
#endif
#if defined(CAHNHILLIARD)
                                     , multiphase_ &
#endif
                                     )
         USE HexMeshClass
         use FluidData
         IMPLICIT NONE
         CLASS(HexMesh)                      :: mesh
#if defined(NAVIERSTOKES) || defined(INCNS)
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
         type(Multiphase_t),     intent(in)  :: multiphase_
#endif
      END SUBROUTINE UserDefinedFinalSetup_f

      subroutine UserDefinedInitialCondition_f(mesh &
#if defined(NAVIERSTOKES) || defined(INCNS)
                                     , thermodynamics_ &
                                     , dimensionless_  &
                                     , refValues_ & 
#endif
#if defined(CAHNHILLIARD)
                                     , multiphase_ &
#endif
                                     )
         use smconstants
         use physicsstorage
         use hexmeshclass
         use fluiddata
         implicit none
         class(hexmesh)                      :: mesh
#if defined(NAVIERSTOKES) || defined(INCNS)
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
         type(Multiphase_t),     intent(in)  :: multiphase_
#endif
      end subroutine UserDefinedInitialCondition_f
#if defined(NAVIERSTOKES) || defined(INCNS)
      subroutine UserDefinedState_f(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
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
         real(kind=RP), intent(inout)  :: Q(NNS)
         type(Thermodynamics_t),    intent(in)  :: thermodynamics_
         type(Dimensionless_t),     intent(in)  :: dimensionless_
         type(RefValues_t),         intent(in)  :: refValues_
      end subroutine UserDefinedState_f

      subroutine UserDefinedNeumann_f(x, t, nHat, U_x, U_y, U_z)
         use SMConstants
         use PhysicsStorage
         use FluidData
         implicit none
         real(kind=RP), intent(in)     :: x(NDIM)
         real(kind=RP), intent(in)     :: t
         real(kind=RP), intent(in)     :: nHat(NDIM)
         real(kind=RP), intent(inout)  :: U_x(NGRADNS)
         real(kind=RP), intent(inout)  :: U_y(NGRADNS)
         real(kind=RP), intent(inout)  :: U_z(NGRADNS)
      end subroutine UserDefinedNeumann_f
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedPeriodicOperation_f(mesh, time, dt, Monitors)
         use SMConstants
         USE HexMeshClass
         use MonitorsClass
         IMPLICIT NONE
         CLASS(HexMesh)               :: mesh
         REAL(KIND=RP)                :: time
         REAL(KIND=RP)                :: dt
         type(Monitor_t), intent(in) :: monitors
      END SUBROUTINE UserDefinedPeriodicOperation_f
!
!//////////////////////////////////////////////////////////////////////// 
! 
#if defined(NAVIERSTOKES) || defined(INCNS)
      subroutine UserDefinedSourceTermNS_f(x, Q, time, S, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         USE HexMeshClass
         use FluidData
         use PhysicsStorage
         IMPLICIT NONE
         real(kind=RP),             intent(in)  :: x(NDIM)
         real(kind=RP),             intent(in)  :: Q(NNS)
         real(kind=RP),             intent(in)  :: time
         real(kind=RP),             intent(out) :: S(NNS)
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
      end subroutine UserDefinedSourceTermNS_f
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedFinalize_f(mesh, time, iter, maxResidual &
#if defined(NAVIERSTOKES) || defined(INCNS)
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
         use SMConstants
         USE HexMeshClass
         use FluidData
         use MonitorsClass
         IMPLICIT NONE
         CLASS(HexMesh)                        :: mesh
         REAL(KIND=RP)                         :: time
         integer                               :: iter
         real(kind=RP)                         :: maxResidual
#if defined(NAVIERSTOKES) || defined(INCNS)
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
      END SUBROUTINE UserDefinedFinalize_f

      SUBROUTINE UserDefinedTermination_f
         implicit none
      END SUBROUTINE UserDefinedTermination_f
   end interface
   
end module ProblemFileFunctions

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
#if defined(NAVIERSTOKES) || defined(INCNS)
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
            CLASS(HexMesh)                      :: mesh
#if defined(NAVIERSTOKES) || defined(INCNS)
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
#if defined(NAVIERSTOKES) || defined(INCNS)
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
#if defined(NAVIERSTOKES) || defined(INCNS)
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
            real(kind=RP)  :: qq, u, v, w, p, x(NDIM), vfcn
#if defined(NAVIERSTOKES)
            real(kind=RP)  :: Q(NNS), phi, theta
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
      
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 
                  qq = 1.0_RP
                  u  = qq*cos(theta)*cos(phi)
                  v  = qq*sin(theta)*cos(phi)
                  w  = qq*sin(phi)
      
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
!           ------------------------------------------------------
!           Incompressible Navier-Stokes default initial condition
!           ------------------------------------------------------
!
#if defined(INCNS)
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 
                  x = mesh % elements(eID) % geom % X(:,i,j,k)
                  vfcn = cos(PI*(sum(x)))
                  mesh % elements(eID) % storage % q(:,i,j,k) = [1.0_RP,vfcn,-2.0_RP*vfcn,vfcn, 2.0_RP*vfcn - 3.0_RP*dimensionless_ % mu(1)*PI*sin(PI*sum(x))] 
               end do;        end do;        end do
               end associate
            end do
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
#if defined(NAVIERSTOKES) || defined(INCNS)
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
            real(kind=RP), intent(inout)  :: Q(NNS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
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
!           ----------------------------------------------------------
!           Called before every time-step to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use MonitorsClass
            IMPLICIT NONE
            CLASS(HexMesh)               :: mesh
            REAL(KIND=RP)                :: time
            REAL(KIND=RP)                :: dt
            type(Monitor_t), intent(in) :: monitors
            
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
#if defined(NAVIERSTOKES) || defined(INCNS)
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
            real(kind=RP),             intent(in)  :: Q(NNS)
            real(kind=RP),             intent(in)  :: time
            real(kind=RP),             intent(out) :: S(NNS)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
!
!           ---------------
!           Local variables
!           ---------------
!
            integer  :: i, j, k, eID
            real(kind=RP) :: cos_, sin_

#ifdef INCNS

            cos_ = cos(PI*(sum(x)-2.0_RP*time))
            sin_ = sin(PI*(sum(x)-2.0_RP*time))
            
            S(INSRHO) = 0.0_RP
            S(INSRHOU) = 0.0_RP
            S(INSRHOV) = -3.0_RP*PI*(2.0_RP*sin_ + 3.0_RP*dimensionless_ % mu(1)*PI*cos_)
            S(INSRHOW) = 0.0_RP
            S(INSP)    = 2.0_RP*PI*(2.0_RP*sin_ + 3.0_RP*PI*dimensionless_ % mu(1)*cos_)

!           S(:) = x(1) + x(2) + x(3) + time
#endif
   
         end subroutine UserDefinedSourceTermNS
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual &
#if defined(NAVIERSTOKES) || defined(INCNS)
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
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
#if defined(NAVIERSTOKES) || defined(INCNS)
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
            integer,       parameter  :: no_of_iters = 331
            real(kind=RP), parameter  :: saved_errors(5) = [ 1.123289361943051E-005_RP, &
                                                             1.400770155204512E-003_RP, &
                                                             1.979726335452486E-003_RP, &
                                                             1.400770155204491E-003_RP, &
                                                             4.479970674933102E-003_RP]
            real(kind=RP), parameter  :: saved_residuals(5) = [1.107385575522812E-003_RP, &
                                                               6.29101636426761_RP, &
                                                               12.5774994950143_RP, &
                                                               6.29101636426761_RP, &
                                                               12.6895468014350_RP]
            CHARACTER(LEN=29)                  :: testName = "Incompressible convergence"
            real(kind=RP)                      :: error(5)
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
            integer                            :: eID, i, j, k, eq, fid,N
            real(kind=RP)                      :: locErr(NCONS), L2Err(NCONS), Q(NCONS), x(NDIM), rho
            real(kind=RP)                      :: wi, wj, wk, Jijk
            real(kind=RP)                      :: cos_, sin_
            real(kind=RP), parameter           :: w_LGL(0:5) = [0.066666666666667_RP, &
                                                      0.378474956297847_RP, &
                                                      0.554858377035486_RP, &
                                                      0.554858377035486_RP, &
                                                      0.378474956297847_RP, &
                                                      0.066666666666667_RP ]
#ifdef INCNS
!
!           *********************************
!           Check the L-inf norm of the error
!           *********************************
!
            L2Err = 0.0_RP

            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 
                  x = mesh % elements(eID) % geom % X(:,i,j,k)

                  cos_ = cos(PI*(sum(x)-2.0_RP*time))
                  sin_ = sin(PI*(sum(x)-2.0_RP*time))
                  Q = mesh % elements(eID) % storage % q(:,i,j,k) 
                  locErr = (Q-[1.0_RP,cos_,-2.0_RP*cos_,cos_,2.0_RP*cos_-3.0_RP*dimensionless_ % mu(1)*PI*sin_])**2

                  wi = w_LGL(i)
                  wj = w_LGL(j)
                  wk = w_LGL(k)
                  Jijk = mesh % elements(eID) % geom % jacobian(i,j,k)

                  L2Err = L2Err + wi*wj*wk*Jijk*locErr 

               end do;        end do;        end do
               N = Nz
               end associate
            end do

            L2Err = sqrt(L2Err)

            print*, "ERROR = ", L2Err


            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()

            CALL FTAssertEqual(expectedValue =  no_of_iters, &
                               actualValue   =  iter, &
                               msg           = "Number of time steps")

            
            CALL FTAssertEqual(expectedValue = saved_errors(1)+1.0_RP, &
                               actualValue   = L2Err(1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Density error")

            CALL FTAssertEqual(expectedValue = saved_errors(2)+1.0_RP, &
                               actualValue   = L2Err(2)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "X-Momentum error")

            CALL FTAssertEqual(expectedValue = saved_errors(3)+1.0_RP, &
                               actualValue   = L2Err(3)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Y-Momentum error")

            CALL FTAssertEqual(expectedValue = saved_errors(4)+1.0_RP, &
                               actualValue   = L2Err(4)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Z-Momentum error")

            CALL FTAssertEqual(expectedValue = saved_errors(5)+1.0_RP, &
                               actualValue   = L2Err(5)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Pressure error")

            CALL FTAssertEqual(expectedValue = saved_residuals(1)+1.0_RP, &
                               actualValue   = monitors % residuals % values(1,1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Density residual")

            CALL FTAssertEqual(expectedValue = saved_residuals(2)+1.0_RP, &
                               actualValue   = monitors % residuals % values(2,1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "X-Momentum residual")

            CALL FTAssertEqual(expectedValue = saved_residuals(3)+1.0_RP, &
                               actualValue   = monitors % residuals % values(3,1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Y-Momentum residual")

            CALL FTAssertEqual(expectedValue = saved_residuals(4)+1.0_RP, &
                               actualValue   = monitors % residuals % values(4,1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Z-Momentum residual")

            CALL FTAssertEqual(expectedValue = saved_residuals(5)+1.0_RP, &
                               actualValue   = monitors % residuals % values(5,1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Pressure residual")


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
      