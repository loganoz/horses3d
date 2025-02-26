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
            real(kind=RP) :: xHalf,yHalf,zHalf,&
                              tempPara01,tempPara02,tempPara03, &
                              gaussionWave
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
#if defined(SCALAR)
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz
                  do j = 0, Ny
                     do i = 0, Nx 
                        x = mesh % elements(eID) % geom % x(:,i,j,k)
                        gaussionWave = sin(x(IX)*PI*2)
                        ! gaussionWave = exp(x(IX)*PI*2)
                        ! mesh % elements(eID) % storage % q(:,i,j,k) = [x(IX)]
                        ! mesh % elements(eID) % storage % q(:,i,j,k) = [1_RP*gaussionWave+10_RP]
                        ! mesh % elements(eID) % storage % q(:,i,j,k) = [3.0_RP]
                        mesh % elements(eID) % storage % q(:,i,j,k) = -[0*1_RP*gaussionWave+100_RP]
                        ! mesh % elements(eID) % storage % q(:,i,j,k) = [10_RP*exp(gaussionWave)] 
                        ! mesh % elements(eID) % storage % q(:,i,j,k) = [1.0_RP] 
                     end do
                  end do
               end do
               end associate
            end do
#endif
#if defined(SCALAR_INS_V04)
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz
                  do j = 0, Ny
                     do i = 0, Nx 
                        x = mesh % elements(eID) % geom % x(:,i,j,k)
                        gaussionWave = sin(x(IX)*PI*2)
                        ! gaussionWave = exp(x(IX)*PI*2)
                        ! mesh % elements(eID) % storage % q(:,i,j,k) = [x(IX)]
                        ! mesh % elements(eID) % storage % q(:,i,j,k) = [1_RP*gaussionWave+10_RP]
                        ! mesh % elements(eID) % storage % q(:,i,j,k) = [3.0_RP]
                        mesh % elements(eID) % storage % q(:,i,j,k) = -1_RP*gaussionWave+100_RP
                        ! mesh % elements(eID) % storage % q(:,i,j,k) = [10_RP*exp(gaussionWave)] 
                        ! mesh % elements(eID) % storage % q(:,i,j,k) = [1.0_RP] 
                     end do
                  end do
               end do
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
#if defined(NAVIERSTOKES) || defined(INCNS) || defined(SCALAR) || defined(SCALAR_INS_V04)
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
            real(kind=RP),             intent(out) :: S(NCONS)
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

            ! x = mesh % elements(eID) % geom % X(:,i,j,k)
                           
            ! cos_ = -300_RP
            ! cos_ = -100*sin(2.5_RP*PI*x(1))
            ! cos_ = -12/( cos(2*x(1)) * cos(2*x(1)) )
            ! cos_ = -2_RP * PI * PI * sin(  PI*x(1) ) * sin( PI*x(3) )
            cos_ = -8_RP * PI * PI * sin(2.0_RP * PI*x(1) ) * sin(2.0_RP * PI*x(3) )
            ! cos_ = -2_RP * PI * PI * sin(1.0_RP* PI*x(1) ) * sin(1_RP* PI*x(3) )
            ! cos_ = -12_RP * PI * PI * sin( 2*PI*x(1) ) * sin( 2*PI*x(2) ) * sin( 2*PI*x(3) )
            ! cos_ = -8_RP * PI * PI * sin( 2*PI*x(1) )   * sin( 2*PI*x(3) )
            do i = 1, NCONS  
               S(i) = cos_
               ! S(i) = 2._RP
            enddo
                           
   
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
            use NodalStorageClass

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
            real(kind=RP)                      :: locErr(NCONS), L2Err(NCONS), L2ErrRelative(NCONS), Q(NCONS), x(NDIM), rho
            real(kind=RP)                      :: wi, wj, wk, Jijk
            real(kind=RP)                      :: cos_, sin_
            type(NodalStorage_t), target  :: spA(3)       ! Nodal storage in the different directions for e_plus  - local copy
            type(NodalStorage_t), pointer :: sp1
            type(NodalStorage_t), pointer :: sp2
            type(NodalStorage_t), pointer :: sp3 
            ! real(kind=RP), parameter           :: w_LGL(0:1) = [1.0_RP,1.0_RP] 
! 
            ! real(kind=RP), parameter           :: w_LGL(0:2) = [0.33333333333333331, &
            !                                                    1.3333333333333333,   &
            !                                                    0.33333333333333331]

            ! real(kind=RP), parameter           :: w_LGL(0:3) = [  0.16666666666666666, &
            !                                                       0.83333333333333337, &
            !                                                       0.83333333333333337, &
            !                                                       0.16666666666666666 &
            !                                                    ]
            real(kind=RP), parameter           :: w_LGL(0:4) = [  0.10000000000000001, &       
                                                                  0.54444444444444440, &       
                                                                  0.71111111111111114, &       
                                                                  0.54444444444444440, &       
                                                                  0.10000000000000001  &
                                                               ]     
                                                               
            ! real(kind=RP), parameter           :: w_LGL(0:5) = [  0.066666666666667_RP, &
            !                                                       0.378474956297847_RP, &
            !                                                       0.554858377035486_RP, &
            !                                                       0.554858377035486_RP, &
            !                                                       0.378474956297847_RP, &
            !                                                       0.066666666666667_RP &
            !                                                       ]

            ! 0.0666666666666667 0.3784749562978469 0.5548583770354865
!  0.5548583770354865 0.3784749562978469 0.0666666666666667
            ! 0.06666667 0.37847496 0.55485838 0.55485838 0.37847496 0.06666667

            ! real(kind=RP), parameter           ::  w_LGL(0:6) = [  4.7619047619047616E-002, &  
            !                                                       0.27682604736156574    , &   
            !                                                       0.43174538120986261    , &   
            !                                                       0.48761904761904762    , &   
            !                                                       0.43174538120986261    , &   
            !                                                       0.27682604736156574    , &    
            !                                                       4.7619047619047616E-002  &
            !                                                       ] 

            ! real(kind=RP), parameter           :: w_LGL(0:7) = [  3.5714285714285712E-002 , &  
            !                                                       0.21070422714350615     , &   
            !                                                       0.34112269248350441     , &   
            !                                                       0.41245879465870372     , &   
            !                                                       0.41245879465870372     , &   
            !                                                       0.34112269248350441     , &   
            !                                                       0.21070422714350615     , &    
            !                                                       3.5714285714285712E-002   &
            !                                                    ]  
            ! real(kind=RP), parameter           :: w_LGL(0:8) = [  2.7777777777777776E-002, &  
            !                                                       0.16549536156080558    , &   
            !                                                       0.27453871250016165    , &   
            !                                                       0.34642851097304617    , &   
            !                                                       0.37151927437641724    , &   
            !                                                       0.34642851097304617    , &   
            !                                                       0.27453871250016165    , &   
            !                                                       0.16549536156080558    , &    
            !                                                       2.7777777777777776E-002 &
            !                                                    ]  
            ! real(kind=RP), parameter           :: w_LGL(0:9) = [  2.2222222222222223E-002, &  
            !                                                       0.13330599085107009    , &   
            !                                                       0.22488934206312652    , &   
            !                                                       0.29204268367968378    , &   
            !                                                       0.32753976118389744    , &   
            !                                                       0.32753976118389744    , &   
            !                                                       0.29204268367968378    , &   
            !                                                       0.22488934206312652    , &   
            !                                                       0.13330599085107009    , &    
            !                                                       2.2222222222222223E-002 &
            !                                                    ]
            ! real(kind=RP), parameter           :: w_LGL(0:10) = [ 1.8181818181818181E-002, & 
            !                                                       0.10961227326699503    , &  
            !                                                       0.18716988178030519    , &  
            !                                                       0.24804810426402829    , &  
            !                                                       0.28687912477900784    , &  
            !                                                       0.30021759545569071    , &  
            !                                                       0.28687912477900784    , &  
            !                                                       0.24804810426402829    , &  
            !                                                       0.18716988178030519    , &  
            !                                                       0.10961227326699503    , &   
            !                                                       1.8181818181818181E-002  &
            !                                                    ]
         ! real(kind=RP), parameter           :: w_LGL(0:15) = [  & 
         !                         8.3333333333333332E-003, &  
         !                         5.0850361005920025E-002, &  
         !                         8.9393697325930860E-002, & 
         !                         0.12425538213251405    , &  
         !                         0.15402698080716426    , &  
         !                         0.17749191339170411    , &  
         !                         0.19369002382520362    , &  
         !                         0.20195830817823002    , &  
         !                         0.20195830817823002    , &  
         !                         0.19369002382520362    , &  
         !                         0.17749191339170411    , &  
         !                         0.15402698080716426    , &  
         !                         0.12425538213251405    , &   
         !                         8.9393697325930860E-002, &  
         !                         5.0850361005920025E-002, &  
         !                         8.3333333333333332E-003 &
         !                                  ]
            ! real(kind=RP), parameter           :: w_LGL(0:20) = [  & 
            !                                                       4.7619047619047623E-003,&   
            !                                                       2.9184840098505565E-002,&   
            !                                                       5.1843169000849641E-002,&   
            !                                                       7.3273918185074172E-002,&   
            !                                                       9.2985467957885995E-002,&  
            !                                                       0.11051708321912326    ,&   
            !                                                       0.12545812119086891    ,&   
            !                                                       0.13745846286004129    ,&   
            !                                                       0.14623686244797737    ,&   
            !                                                       0.15158757511168136    ,&   
            !                                                       0.15338519033217496    ,&   
            !                                                       0.15158757511168136    ,&   
            !                                                       0.14623686244797737    ,&   
            !                                                       0.13745846286004129    ,&   
            !                                                       0.12545812119086891    ,&   
            !                                                       0.11051708321912326    ,&    
            !                                                       9.2985467957885995E-002,&   
            !                                                       7.3273918185074172E-002,&   
            !                                                       5.1843169000849641E-002,&   
            !                                                       2.9184840098505565E-002,&   
            !                                                       4.7619047619047623E-003 &
            !                                                    ]
! ==========================================================================================================
            ! real(kind=RP), parameter           :: w_LGL(0:1) = [  & 
            !                                                    1.0000000000000000 , &       
            !                                                    1.0000000000000000   &
            !                                                    ]
            ! real(kind=RP), parameter           :: w_LGL(0:2) = [  & 
            ! 0.55555555555555569     ,&  
            ! 0.88888888888888884     ,&  
            ! 0.55555555555555569     &
            ! ]
            
            ! real(kind=RP), parameter           :: w_LGL(0:3) = [  & 
            ! 0.34785484513745385, &       
            ! 0.65214515486254632, &       
            ! 0.65214515486254632, &       
            ! 0.34785484513745385 &     
            ! ]

            ! real(kind=RP), parameter           :: w_LGL(0:4) = [  & 
            !                                         0.23692688505618911    , &   
            !                                         0.47862867049936669    , &   
            !                                         0.56888888888888889    , &   
            !                                         0.47862867049936669    , &   
            !                                         0.23692688505618911     & 
            !                                        ]
            ! real(kind=RP), parameter           :: w_LGL(0:5) = [  & 
            !                                              0.17132449237917019     ,&  
            !                                              0.36076157304813833     ,&  
            !                                              0.46791393457269093     ,&  
            !                                              0.46791393457269093     ,&  
            !                                              0.36076157304813833     ,&  
            !                                              0.17132449237917019     &
            !                                              ]
            ! real(kind=RP), parameter           :: w_LGL(0:6) = [  & 
            !                                     0.12948496616886979  ,&     
            !                                     0.27970539148927670  ,&     
            !                                     0.38183005050511903  ,&     
            !                                     0.41795918367346940  ,&     
            !                                     0.38183005050511903  ,&     
            !                                     0.27970539148927670  ,&     
            !                                     0.12948496616886979 &
            !                                     ]     

            ! real(kind=RP), parameter           :: w_LGL(0:7) = [  &
            !                                              0.10122853629037630  , &     
            !                                              0.22238103445337470  , &     
            !                                              0.31370664587788727  , &     
            !                                              0.36268378337836193  , &     
            !                                              0.36268378337836193  , &     
            !                                              0.31370664587788727  , &     
            !                                              0.22238103445337470  , &     
            !                                              0.10122853629037630   &
            !                                              ]   
            ! real(kind=RP), parameter           :: w_LGL(0:8) = [  & 
            !                                                 8.1274388361574565E-002, & 
            !                                                 0.18064816069485751    , &  
            !                                                 0.26061069640293549    , &  
            !                                                 0.31234707704000264    , &  
            !                                                 0.33023935500125978    , &  
            !                                                 0.31234707704000264    , &  
            !                                                 0.26061069640293549    , &  
            !                                                 0.18064816069485751    , &   
            !                                                 8.1274388361574565E-002 &
            !                                                 ]
         ! real(kind=RP), parameter           :: w_LGL(0:8) = [  &  
         !                                                 8.1274388361574565E-002 ,  &
         !                                                 0.18064816069485751     ,  & 
         !                                                 0.26061069640293549     ,  & 
         !                                                 0.31234707704000264     ,  & 
         !                                                 0.33023935500125978     ,  & 
         !                                                 0.31234707704000264     ,  & 
         !                                                 0.26061069640293549     ,  & 
         !                                                 0.18064816069485751     ,  &  
         !                                                 8.1274388361574565E-002  &
         !                                                 ]
         ! real(kind=RP), parameter           :: w_LGL(0:9) = [  &  
         !                                              6.6671344308688527E-002  ,&
         !                                              0.14945134915058061      ,& 
         !                                              0.21908636251598199      ,& 
         !                                              0.26926671930999607      ,& 
         !                                              0.29552422471475270      ,& 
         !                                              0.29552422471475270      ,& 
         !                                              0.26926671930999607      ,& 
         !                                              0.21908636251598199      ,& 
         !                                              0.14945134915058061      ,&  
         !                                              6.6671344308688527E-002 &
         !                                              ]
                                                            
         ! real(kind=RP), parameter           :: w_LGL(0:10) = [  & 
         !                                                        5.5668567116173601E-002,&
         !                                                        0.12558036946490472,&
         !                                                        0.18629021092773429 ,&     
         !                                                        0.23319376459199040     ,&       
         !                                                        0.26280454451024682     ,&       
         !                                                        0.27292508677790062     ,&       
         !                                                        0.26280454451024682     ,&       
         !                                                        0.23319376459199040     ,&       
         !                                                        0.18629021092773429     ,&       
         !                                                        0.12558036946490472     ,&        
         !                                                        5.5668567116173601E-002 &
         !                                                        ]
         ! real(kind=RP), parameter           :: w_LGL(0:20) = [  &   
         !                                                  1.6017228257774130E-002, &  
         !                                                  3.6953789770852459E-002, &  
         !                                                  5.7134425426857260E-002, &  
         !                                                  7.6100113628379415E-002, &  
         !                                                  9.3444423456033807E-002, & 
         !                                                  0.10879729916714852    , &  
         !                                                  0.12183141605372855    , &  
         !                                                  0.13226893863333730    , &  
         !                                                  0.13988739479107301    , &  
         !                                                  0.14452440398996996    , &  
         !                                                  0.14608113364969041    , &  
         !                                                  0.14452440398996996    , &  
         !                                                  0.13988739479107301    , &  
         !                                                  0.13226893863333730    , &  
         !                                                  0.12183141605372855    , &  
         !                                                  0.10879729916714852    , &   
         !                                                  9.3444423456033807E-002, &  
         !                                                  7.6100113628379415E-002, &  
         !                                                  5.7134425426857260E-002, &  
         !                                                  3.6953789770852459E-002, &  
         !                                                  1.6017228257774130E-002  &
         !                                                  ]
            ! type(NodalStorage_t)          :: spA(3)       ! Nodal storage in the different directions for e_plus  - local copy
 

#ifdef SCALAR
!
!           *********************************
!           Check the L-inf norm of the error
!           *********************************
!
            L2Err = 0.0_RP
            L2ErrRelative = 0.0_RP

            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )


               do k = 0, Nz;  
                  do j = 0, Ny;  
                     do i = 0, Nx 
                        x = mesh % elements(eID) % geom % X(:,i,j,k)

                        cos_ = -sin( 2*PI*x(1) ) * sin( 2*PI*x(3) ) + 0.0_RP
                        ! sin_ = sin(PI*(sum(x)-2.0_RP*time))
                        Q = mesh % elements(eID) % storage % q(:,i,j,k) 
                        locErr = (Q-[cos_])*(Q-[cos_])


                        wi = w_LGL(i)
                        ! wj = w_LGL(j)
                        wj = 1
                        wk = w_LGL(k)
                        Jijk = mesh % elements(eID) % geom % jacobian(i,j,k)
                        write (*,*) "Q, cos_, Q-cos_,Jijk =============", eID,i,j,k,Q, cos_, Q-cos_,Jijk
                        write (*,*) "wi*wj*wk,locErr,L2Err =============", wi*wj*wk,locErr,L2Err

                        L2Err = L2Err + wi*wj*wk*Jijk*locErr
                        L2ErrRelative = L2ErrRelative +  wi*wj*wk*Jijk* [cos_*cos_]
                        ! L2ErrRelative = L2ErrRelative +  wi*wj*wk*Jijk*locErr/([cos_*cos_])
                        ! L2Err = L2Err + sp1%w(i)*sp2%w(j)*sp3%w(k)*Jijk*locErr 

                     end do;        
                  end do;        
               end do
               N = Nz
               end associate
            end do

            L2Err = sqrt(L2Err)

            L2ErrRelative = L2Err / sqrt(L2ErrRelative)


            print*, "ERROR = ", L2Err
            print*, "Relative error = ",  L2ErrRelative

            write(*,*) "w_LGL,size(w_LGL), p, absolute, relative error  = ",  &
                                                            w_LGL(0) , size(w_LGL)   ,               &
                                                            mesh%no_of_elements,          &
                                                            (mesh%no_of_elements)**(1._RP/3),          &
                                                            mesh%faces(5)%geom%h,         &
                                                            mesh%NDOF,                    &
                                                            mesh % elements(1) % Nxyz(1), &
                                                            L2Err,                        &
                                                            L2ErrRelative






#endif

#ifdef SCALAR_INS_V04
!
!           *********************************
!           Check the L-inf norm of the error
!           *********************************
!
            L2Err = 0.0_RP
            L2ErrRelative = 0.0_RP

            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )


               do k = 0, Nz;  
                  do j = 0, Ny;  
                     do i = 0, Nx 
                        x = mesh % elements(eID) % geom % X(:,i,j,k)

                        cos_ = -sin( 2*PI*x(1) ) * sin( 2*PI*x(3) ) + 0.0_RP
                        ! sin_ = sin(PI*(sum(x)-2.0_RP*time))
                        Q = mesh % elements(eID) % storage % q(:,i,j,k) 
                        locErr = (Q-[cos_,cos_,cos_,cos_])*(Q-[cos_,cos_,cos_,cos_])


                        wi = w_LGL(i)
                        ! wj = w_LGL(j)
                        wj = 1
                        wk = w_LGL(k)
                        Jijk = mesh % elements(eID) % geom % jacobian(i,j,k)
                        write (*,*) "Q, cos_, Q-cos_,Jijk =============", eID,i,j,k,Q, cos_, Q-cos_,Jijk
                        write (*,*) "wi*wj*wk,locErr,L2Err =============", wi*wj*wk,locErr,L2Err

                        L2Err = L2Err + wi*wj*wk*Jijk*locErr
                        L2ErrRelative = L2ErrRelative +  wi*wj*wk*Jijk* [cos_*cos_,cos_*cos_,cos_*cos_,cos_*cos_]
                        ! L2ErrRelative = L2ErrRelative +  wi*wj*wk*Jijk*locErr/([cos_*cos_])
                        ! L2Err = L2Err + sp1%w(i)*sp2%w(j)*sp3%w(k)*Jijk*locErr 

                     end do;        
                  end do;        
               end do
               N = Nz
               end associate
            end do

            L2Err = sqrt(L2Err)

            L2ErrRelative = L2Err / sqrt(L2ErrRelative)


            print*, "ERROR = ", L2Err
            print*, "Relative error = ",  L2ErrRelative

            write(*,*) "w_LGL,size(w_LGL), p, absolute, relative error  = ",  &
                                                            w_LGL(0) , size(w_LGL)   ,               &
                                                            mesh%no_of_elements,          &
                                                            (mesh%no_of_elements)**(1._RP/3),          &
                                                            mesh%faces(5)%geom%h,         &
                                                            mesh%NDOF,                    &
                                                            mesh % elements(1) % Nxyz(1), &
                                                            L2Err,                        &
                                                            L2ErrRelative






#endif

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
               do k = 0, Nz;  
                  do j = 0, Ny;  
                     do i = 0, Nx 
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

                        end do;        
                     end do;        
                  end do
                  N = Nz
               end associate
            end do

            L2Err = sqrt(L2Err)

            print*, "ERROR = ", L2Err





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
      
