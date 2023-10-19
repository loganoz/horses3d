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

#if defined(NAVIERSTOKES) & !(defined(SPALARTALMARAS))
module Convergence_analysis
   use SMConstants
   use PhysicsStorage
   implicit none

   private 
   public   state_source_in_point
   contains
    subroutine state_source_in_point(x,y,z,t,gm1,gMa2,mu,kap,Q,S, QDot)
!
!      ********************************************
!      Implementation of a source term a la Diego
!      ********************************************
!
       implicit none
       real(kind=RP), intent(in)  :: x,y,z,t,gm1,gMa2,mu,kap
       real(kind=RP), intent(out) :: Q(NCONS), S(NCONS), QDot(NCONS)
!
!      ---------------
!      Local variables
!      ---------------
!
       real(kind=RP), parameter :: c_u = 1.0_RP, c_v = 2.0_RP, c_w = -c_u-c_v
       real(kind=RP), parameter :: S_div_TRef_Sutherland = 0.381923076923077_RP
       real(kind=RP) :: arg, dargdx, dargdt, rho, u, v, w, e, p, Temp, vtot
       real(kind=RP) :: invRho
       real(kind=RP) :: rho_x(NDIM), u_x(NDIM,1), v_x(NDIM,1), w_x(NDIM,1), e_x(NDIM), p_x(NDIM), Temp_x(NDIM)
       real(kind=RP) :: rho_xx(NDIM,NDIM), u_xx(NDIM,NDIM), v_xx(NDIM,NDIM), w_xx(NDIM,NDIM), e_xx(NDIM,NDIM), Temp_xx(NDIM,NDIM)
       real(kind=RP) :: tau(NDIM,NDIM), grad_u(NDIM,NDIM), eye3(NDIM,NDIM), div_tau(NDIM)
       real(kind=RP) :: rho_t, u_t, v_t, w_t, e_t, p_t
       real(kind=RP) :: div_v, suther, dsutherdT, visc_work, heat_flux

       eye3 = 0.0_RP
       eye3(1,1) = 1.0_RP
       eye3(2,2) = 1.0_RP
       eye3(3,3) = 1.0_RP
        
       arg = PI*(x+y+z-2.0_RP*t)
       dargdx = PI
       dargdt = -2.0_RP*PI
!
!      State 
!      -----
       rho       = 2.0_RP + 0.1_RP*sin(arg)
       invRho    = 1.0_RP / rho
       u         = c_u*rho - 1.0_RP
       v         = c_v*rho - 4.0_RP
       w         = c_w*rho + 5.0_RP
       e         = rho
       vtot      = u*u+v*v+w*w
       p         = gm1*(rho*e-0.5_RP*rho*vtot)
       Temp      = gMa2*gm1*(e - 0.5_RP*vtot)
       suther    = (1._RP + S_div_TRef_Sutherland)/(Temp + S_div_TRef_Sutherland)*Temp*SQRT(Temp)
       dsutherdT = 1.5_RP*suther/Temp - suther/(S_div_TRef_Sutherland+Temp)

       Q = [rho, rho*u, rho*v, rho*w, rho*e]
!
!      Time derivatives
!      ----------------
       rho_t = 0.1_RP*dargdt*cos(arg)
       u_t   = c_u*rho_t
       v_t   = c_v*rho_t
       w_t   = c_w*rho_t
       e_t   = rho_t

       QDot = [rho_t,rho_t*u + rho*u_t, rho_t*v + rho*v_t, rho_t*w + rho*w_t, rho_t*e + rho*e_t]
!
!      Gradients
!      ---------
       rho_x = 0.1_RP*dargdx*cos(arg)
       u_x(:,1)   = c_u*rho_x
       v_x(:,1)   = c_v*rho_x
       w_x(:,1)   = c_w*rho_x
       e_x   = rho_x
       p_x   = gm1*rho_x*(e-0.5_RP*vtot) + gm1*rho*(e_x-u*u_x(:,1)-v*v_x(:,1)-w*w_x(:,1))
       Temp_x   = gMa2*gm1*(e_x - u*u_x(:,1) - v*v_x(:,1) - w*w_x(:,1))
       div_v = u_x(IX,1) + v_x(IY,1) + w_x(IZ,1)

       grad_u(:,1) = u_x(:,1)
       grad_u(:,2) = v_x(:,1)
       grad_u(:,3) = w_x(:,1)

       tau = mu*suther*(grad_u + transpose(grad_u) - 2.0_RP*div_v*eye3/3.0_RP)
!
!      Hessians
!      --------
       rho_xx = -0.1_RP*POW2(dargdx)*sin(arg)
       u_xx   = c_u*rho_xx
       v_xx   = c_v*rho_xx
       w_xx   = c_w*rho_xx
       e_xx   = rho_xx
       Temp_xx   = gMa2*gm1*(e_xx - matmul(u_x,transpose(u_x))-matmul(v_x,transpose(v_x))-matmul(w_x,transpose(w_x))-u*u_xx-v*v_xx-w*w_xx)

       div_tau = u_xx(1,:) + v_xx(2,:) + w_xx(3,:) + [u_xx(1,1)+u_xx(2,2)+u_xx(3,3), v_xx(1,1)+v_xx(2,2)+v_xx(3,3), w_xx(1,1)+w_xx(2,2)+w_xx(3,3)] &  
                     - (2.0_RP/3.0_RP)*[u_xx(1,1)+v_xx(2,1)+w_xx(3,1),u_xx(2,1)+v_xx(2,2)+w_xx(2,3),u_xx(1,3)+v_xx(2,3)+w_xx(3,3)] 
       div_tau = mu*suther*div_tau


       div_tau = div_tau + (dsutherdT/suther)*matmul(tau,Temp_x)
       visc_work = sum(grad_u*tau) + u*div_tau(IX) + v*div_tau(IY) + w*div_tau(IZ)
       heat_flux = kap*dsutherdT*(POW2(Temp_x(IX)) + POW2(Temp_x(IY)) + POW2(Temp_x(IZ))) + kap*suther*(Temp_xx(IX,IX) + Temp_xx(IY,IY) + Temp_xx(IZ,IZ))
!
!      Source term
!      -----------
       S(IRHO)  = rho_t + u*rho_x(IX) + v*rho_x(IY) + w*rho_x(IZ) + rho*div_v
       S(IRHOU) = u*S(IRHO) + rho*u_t + rho*u*u_x(IX,1) + rho*v*u_x(IY,1) + rho*w*u_x(IZ,1) + p_x(IX) - div_tau(IX)
       S(IRHOV) = v*S(IRHO) + rho*v_t + rho*u*v_x(IX,1) + rho*v*v_x(IY,1) + rho*w*v_x(IZ,1) + p_x(IY) - div_tau(IY)
       S(IRHOW) = w*S(IRHO) + rho*w_t + rho*u*w_x(IX,1) + rho*v*w_x(IY,1) + rho*w*w_x(IZ,1) + p_x(IZ) - div_tau(IZ)
       S(IRHOE) = e*S(IRHO) + rho*e_t + rho*(u*e_x(IX)+v*e_x(IY)+w*e_x(IZ)) + p*div_v + u*p_x(IX)+v*p_x(IY)+w*p_x(IZ) - visc_work - heat_flux

    end subroutine state_source_in_point

end module Convergence_analysis
#endif

module ProblemFileFunctions
   implicit none

   abstract interface
      subroutine UserDefinedStartup_f
      end subroutine UserDefinedStartup_f
   
      SUBROUTINE UserDefinedFinalSetup_f(mesh &
#ifdef FLOW
                                     , thermodynamics_ &
                                     , dimensionless_  &
                                     , refValues_ & 
#endif
#ifdef CAHNHILLIARD
                                     , multiphase_ &
#endif
                                     )
         USE HexMeshClass
         use FluidData
         IMPLICIT NONE
         CLASS(HexMesh)                      :: mesh
#ifdef FLOW
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)  :: multiphase_
#endif
      END SUBROUTINE UserDefinedFinalSetup_f

      subroutine UserDefinedInitialCondition_f(mesh &
#ifdef FLOW
                                     , thermodynamics_ &
                                     , dimensionless_  &
                                     , refValues_ & 
#endif
#ifdef CAHNHILLIARD
                                     , multiphase_ &
#endif
                                     )
         use smconstants
         use physicsstorage
         use hexmeshclass
         use fluiddata
         implicit none
         class(hexmesh)                      :: mesh
#ifdef FLOW
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)  :: multiphase_
#endif
      end subroutine UserDefinedInitialCondition_f
#ifdef FLOW
      subroutine UserDefinedState_f(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
         use SMConstants
         use PhysicsStorage
         use FluidData
         implicit none
         real(kind=RP)  :: x(NDIM)
         real(kind=RP)  :: t
         real(kind=RP)  :: nHat(NDIM)
         real(kind=RP)  :: Q(NCONS)
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
      end subroutine UserDefinedState_f

      subroutine UserDefinedGradVars_f(x, t, nHat, Q, U, GetGradients, thermodynamics_, dimensionless_, refValues_)
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
      end subroutine UserDefinedGradVars_f


      subroutine UserDefinedNeumann_f(x, t, nHat, Q, U_x, U_y, U_z, flux, thermodynamics_, dimensionless_, refValues_)
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
#ifdef FLOW
      subroutine UserDefinedSourceTermNS_f(x, Q, time, S, thermodynamics_, dimensionless_, refValues_ &
#ifdef CAHNHILLIARD
,multiphase_ &
#endif
)
         use SMConstants
         USE HexMeshClass
         use FluidData
         use PhysicsStorage
         IMPLICIT NONE
         real(kind=RP),             intent(in)  :: x(NDIM)
         real(kind=RP),             intent(in)  :: Q(NCONS)
         real(kind=RP),             intent(in)  :: time
         real(kind=RP),             intent(inout) :: S(NCONS)
         type(Thermodynamics_t), intent(in)  :: thermodynamics_
         type(Dimensionless_t),  intent(in)  :: dimensionless_
         type(RefValues_t),      intent(in)  :: refValues_
#ifdef CAHNHILLIARD
         type(Multiphase_t),     intent(in)  :: multiphase_
#endif
      end subroutine UserDefinedSourceTermNS_f
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedFinalize_f(mesh, time, iter, maxResidual &
#ifdef FLOW
                                                 , thermodynamics_ &
                                                 , dimensionless_  &
                                                 , refValues_ & 
#endif   
#ifdef CAHNHILLIARD
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
#ifdef FLOW
         type(Thermodynamics_t), intent(in)    :: thermodynamics_
         type(Dimensionless_t),  intent(in)    :: dimensionless_
         type(RefValues_t),      intent(in)    :: refValues_
#endif
#ifdef CAHNHILLIARD
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
#ifdef FLOW
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#ifdef CAHNHILLIARD
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
#ifdef FLOW
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         subroutine UserDefinedInitialCondition(mesh &
#ifdef FLOW
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#ifdef CAHNHILLIARD
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
#ifdef NAVIERSTOKES
            use Convergence_analysis
#endif
            implicit none
            class(hexmesh)                      :: mesh
#ifdef FLOW
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           local variables
!           ---------------
!
            integer        :: eid, i, j, k
#if defined(NAVIERSTOKES)
            real(kind=RP)  :: Q(NCONS), S(NCONS), phi, theta, x, y, z, QDot(NCONS)
#endif

!
!           ---------------------------------------
!           Navier-Stokes default initial condition
!           ---------------------------------------
!
#if defined(NAVIERSTOKES)
            associate ( gammaM2 => dimensionless_ % gammaM2, &
                        gammaMinus1 => thermodynamics_ % gammaMinus1, &
                        gamma => thermodynamics_ % gamma )
      
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          ny => mesh % elemeNts(eID) % nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 
                  x = mesh % elements(eID) % geom % x(IX,i,j,k)
                  y = mesh % elements(eID) % geom % x(IY,i,j,k)
                  z = mesh % elements(eID) % geom % x(IZ,i,j,k)

                  call state_source_in_point(x,y,z,0.0_RP,gammaMinus1, gammaM2, dimensionless_ % mu, dimensionless_ % kappa, Q,S, QDot)
      
                  mesh % elements(eID) % storage % q(:,i,j,k) = q 
               end do;        end do;        end do
               end associate
            end do

            end associate
#endif
         end subroutine UserDefinedInitialCondition
#ifdef FLOW
         subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
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
#ifdef FLOW
         subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_ &
#ifdef CAHNHILLIARD
, multiphase_ &
#endif
)
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
#ifdef NAVIERSTOKES
            use Convergence_analysis
#endif
            IMPLICIT NONE
            real(kind=RP),             intent(in)  :: x(NDIM)
            real(kind=RP),             intent(in)  :: Q(NCONS)
            real(kind=RP),             intent(in)  :: time
            real(kind=RP),             intent(inout) :: S(NCONS)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           Local variables
!           ---------------
!
            real(kind=RP)  :: Temp(NCONS), Temp2(NCONS)

#ifdef NAVIERSTOKES
!
!           Usage example
!           -------------
            call state_source_in_point( x(IX), x(IY), x(IZ), time, &
                                        thermodynamics_ % gammaMinus1, & 
                                        dimensionless_ % gammaM2, &
                                        dimensionless_ % mu, &
                                        dimensionless_ % kappa, &
                                        Temp, S, Temp2)

#endif
   
         end subroutine UserDefinedSourceTermNS
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual &
#ifdef FLOW
                                                    , thermodynamics_ &
                                                    , dimensionless_  &
                                                    , refValues_ & 
#endif   
#ifdef CAHNHILLIARD
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
#ifdef NAVIERSTOKES
            use Convergence_analysis
#endif
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
#ifdef FLOW
            type(Thermodynamics_t), intent(in)    :: thermodynamics_
            type(Dimensionless_t),  intent(in)    :: dimensionless_
            type(RefValues_t),      intent(in)    :: refValues_
#endif
#ifdef CAHNHILLIARD
            type(Multiphase_t),     intent(in)    :: multiphase_
#endif
            type(Monitor_t),        intent(in)    :: monitors
            real(kind=RP),             intent(in) :: elapsedTime
            real(kind=RP),             intent(in) :: CPUTime
#ifdef NAVIERSTOKES
            real(kind=RP)  :: x, y, z, locErr(5), Q(5), S(5)

            real(kind=RP), parameter  :: residuals(5) = [ 6.2929844029491377E-01_RP, &
                                                          1.8894710750845625E+00_RP, &
                                                          2.5264755519384234E+00_RP, &
                                                          4.4146672499485859E+00_RP, &
                                                          2.5157533917446182E+00_RP]
            real(kind=RP), parameter :: entropyRate = 3.6207469616966779E-07_RP
            real(kind=RP), parameter :: saved_errors(5) = [ 3.3102799903292017E-05_RP, &
                                                            3.6405407003671022E-05_RP, &
                                                            1.5957629488942332E-05_RP, &
                                                            3.4132613560424396E-05_RP, &
                                                            6.6811198249671421E-05_RP]
            real(kind=RP), parameter :: saved_qdot_errors(5) = [ 9.2273153574772720E-04_RP, &
                                                                 9.6259031804949867E-04_RP, &
                                                                 3.8474183309141608E-04_RP, &
                                                                 8.7717012428247660E-04_RP, &
                                                                 1.6750805236722724E-03_RP]
            integer        :: i, j,k, eID 
            CHARACTER(LEN=29)                  :: testName           = "Navier-Stokes convergence"
            real(kind=RP)  :: error(5), QDotError(5), QDot(5), locErrQDot(5)
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
            real(kind=RP), parameter   :: wLGL5(0:5) = [0.066666666666667_RP, &
                                                      0.378474956297847_RP, &
                                                      0.554858377035486_RP, &
                                                      0.554858377035486_RP, &
                                                      0.378474956297847_RP, &
                                                      0.066666666666667_RP ]
            real(kind=RP), parameter  :: wLGL7(0:7) = [  0.035714285714286_RP, &
                                                         0.210704227143506_RP, &
                                                         0.341122692483504_RP, &
                                                         0.412458794658704_RP, &
                                                         0.412458794658704_RP, &
                                                         0.341122692483504_RP, &
                                                         0.210704227143506_RP, &
                                                         0.035714285714286_RP]
!
!           *********************************
!           Check the L-inf norm of the error
!           *********************************
!
            error     = 0.0_RP
            QDotError = 0.0_RP

            do eID = 1, mesh % no_of_elements
               associate(e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  x = e % geom % x(IX,i,j,k)
                  y = e % geom % x(IY,i,j,k)
                  z = e % geom % x(IZ,i,j,k)
 
                  call state_source_in_point(x,y,z,time,thermodynamics_ % gammaMinus1, &
                                             dimensionless_ % mu, dimensionless_ % kappa, dimensionless_ % gammaM2, Q, S, QDot)

                  locErr = e % storage % Q(:,i,j,k) - Q
                  locErrQDot = e % storage % QDot(:,i,j,k) - QDot

                  error = error + e % geom % jacobian(i,j,k)*wLGL7(i)*wLGL7(j)*wLGL7(k)*POW2(locErr)
                  QDotError = QDotError + e % geom % jacobian(i,j,k)*wLGL7(i)*wLGL7(j)*wLGL7(k)*POW2(locErrQDot)

               end do            ; end do ; end do
               end associate
            end do

            error = sqrt(error)
            QDotError = sqrt(QDotError)

!write(STD_OUT,'(A)',advance='no') "["
!do i = 1, 4
!write(STD_OUT,'(ES24.16,A)') monitors % residuals % values(i,1), "_RP, &"
!enddo
!write(STD_OUT,'(ES24.16,A)') monitors % residuals % values(5,1), "_RP]"
!write(STD_OUT,'(A)',advance='no') "["
!do i = 1, 4
!write(STD_OUT,'(ES24.16,A)') error(i), "_RP, &"
!enddo
!write(STD_OUT,'(ES24.16,A)') error(5), "_RP]"
!write(STD_OUT,'(A)',advance='no') "["
!do i = 1, 4
!write(STD_OUT,'(ES24.16,A)') QDotError(i), "_RP, &"
!enddo
!write(STD_OUT,'(ES24.16,A)') QDotError(5), "_RP]"
!
!write(STD_OUT,'(ES24.16,A)') monitors % volumeMonitors(1) % values(1,1), "_RP"


            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()
            
#ifndef _HAS_MPI_
            CALL FTAssertEqual(expectedValue = saved_errors(1)+1.0_RP, &
                               actualValue   = error(1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Continuity error")


            CALL FTAssertEqual(expectedValue = saved_errors(2)+1.0_RP, &
                               actualValue   = error(2)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "X-Momentum error")

            CALL FTAssertEqual(expectedValue = saved_errors(3)+1.0_RP, &
                               actualValue   = error(3)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Y-Momentum error")

            CALL FTAssertEqual(expectedValue = saved_errors(4)+1.0_RP, &
                               actualValue   = error(4)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Z-Momentum error")

            CALL FTAssertEqual(expectedValue = saved_errors(5)+1.0_RP, &
                               actualValue   = error(5)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Energy error")

            CALL FTAssertEqual(expectedValue = saved_qdot_errors(1)+1.0_RP, &
                               actualValue   = QDotError(1)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Continuity qdot error")

            CALL FTAssertEqual(expectedValue = saved_qdot_errors(2)+1.0_RP, &
                               actualValue   = QDotError(2)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "X-Momentum qdot error")

            CALL FTAssertEqual(expectedValue = saved_qdot_errors(3)+1.0_RP, &
                               actualValue   = QDotError(3)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Y-Momentum qdot error")

            CALL FTAssertEqual(expectedValue = saved_qdot_errors(4)+1.0_RP, &
                               actualValue   = QDotError(4)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Z-Momentum qdot error")

            CALL FTAssertEqual(expectedValue = saved_qdot_errors(5)+1.0_RP, &
                               actualValue   = QDotError(5)+1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Energy qdot error")

#endif
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

            CALL FTAssertEqual(expectedValue = entropyRate + 1.0_RP, &
                               actualValue   = monitors % volumeMonitors(1) % values(1,1) + 1.0_RP, &
                               tol           = 1.d-11, &
                               msg           = "Entropy Rate")





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
      