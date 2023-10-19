#include "Includes.h"
module RiemannSolvers_iNSKeywordsModule

     integer,                       parameter :: KEYWORD_LENGTH           = 132
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_SOLVER_NAME_KEY  = "riemann solver"
     character(len=KEYWORD_LENGTH), parameter :: LAMBDA_STABILIZATION_KEY = "lambda stabilization"
     character(len=KEYWORD_LENGTH), parameter :: AVG_NAME_KEY             = "averaging"
!
!    --------------------------
!    Riemann solver definitions
!    --------------------------
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_CENTRAL_NAME = "central"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_LXF_NAME     = "lax-friedrichs"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_EXACT_NAME   = "exact"

     enum, bind(C)
        enumerator :: RIEMANN_CENTRAL = 1, RIEMANN_LXF, RIEMANN_EXACT
     end enum
!
!    -----------------------------
!    Available averaging functions
!    -----------------------------
!
     character(len=KEYWORD_LENGTH), parameter :: STANDARD_AVG_NAME       = "standard"
     character(len=KEYWORD_LENGTH), parameter :: SKEWSYMMETRIC1_AVG_NAME = "skew-symmetric 1"
     character(len=KEYWORD_LENGTH), parameter :: SKEWSYMMETRIC2_AVG_NAME = "skew-symmetric 2"

     enum, bind(C)
        enumerator :: STANDARD_AVG = 1, SKEWSYMMETRIC1_AVG, SKEWSYMMETRIC2_AVG
     end enum

end module RiemannSolvers_iNSKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
!
module RiemannSolvers_iNS
   use SMConstants
   use Physics_iNS
   use PhysicsStorage_iNS
   use VariableConversion_iNS
   use FluidData_iNS

   implicit none

   private
   public whichAverage, whichRiemannSolver
   public SetRiemannSolver, DescribeRiemannSolver
   public RiemannSolver, AveragedStates, TwoPointFlux, ExactRiemannSolver

   abstract interface
      subroutine RiemannSolverFCN(QLeft, QRight, nHat, t1, t2, flux)
         use SMConstants
         use PhysicsStorage_iNS
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
      end subroutine RiemannSolverFCN

      subroutine AveragedStatesFCN(QLeft, QRight, f, g, h)
         use SMConstants
         use PhysicsStorage_iNS
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(out)      :: f(1:NCONS), g(1:NCONS), h(1:NCONS)
      end subroutine AveragedStatesFCN

      subroutine TwoPointFluxFCN(QLeft, QRight, JaL, JaR, fSharp)
         use SMConstants
         use PhysicsStorage_iNS
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
      end subroutine TwoPointFluxFCN
   end interface

   procedure(RiemannSolverFCN),  protected, pointer :: RiemannSolver  => NULL()
   procedure(AveragedStatesFCN), protected, pointer :: AveragedStates => NULL()
   procedure(TwoPointFluxFCN),   protected, pointer :: TwoPointFlux   => NULL()

   integer, protected :: whichRiemannSolver = -1
   integer, protected :: whichAverage = -1
   real(RP)           :: lambdaStab = 1.0_RP
!
!  ========
   contains
!  ========
!
      subroutine SetRiemannSolver(controlVariables)
!
!        -------
!        Modules
!        -------
         use Utilities, only: toLower
         use FTValueDictionaryClass
         use RiemannSolvers_iNSKeywordsModule
!
!        ---------
!        Interface
!        ---------
         type(FTValueDictionary), intent(in) :: controlVariables
!
!        ---------------
!        Local variables
!        ---------------
         character(len=KEYWORD_LENGTH) :: keyword

!
!        --------------------------------------------
!        Choose the Riemann solver (default is exact)
!        --------------------------------------------
         if (controlVariables % containsKey(RIEMANN_SOLVER_NAME_KEY)) then

            keyword = controlVariables % stringValueForKey(RIEMANN_SOLVER_NAME_KEY, KEYWORD_LENGTH)
            call toLower(keyword)

            select case (keyword)
            case(RIEMANN_CENTRAL_NAME)
               RiemannSolver => CentralRiemannSolver
               whichRiemannSolver = RIEMANN_CENTRAL

            case(RIEMANN_LXF_NAME)
               RiemannSolver => LxFRiemannSolver
               whichRiemannSolver = RIEMANN_LXF

            case(RIEMANN_EXACT_NAME)
               RiemannSolver => ExactRiemannSolver
               whichRiemannSolver = RIEMANN_EXACT

            case default
               print*, "Riemann Solver not recognized."
               errorMessage(STD_OUT)
               error stop
            end select

         else
!
!           Select exact by default
!           -----------------------
            RiemannSolver => CentralRiemannSolver
            whichRiemannSolver = RIEMANN_EXACT

         end if
!
!        --------------------
!        Lambda stabilization
!        --------------------
         if (controlVariables % containsKey(LAMBDA_STABILIZATION_KEY)) then
            lambdaStab = controlVariables % doublePrecisionValueForKey(LAMBDA_STABILIZATION_KEY)

         else
!
!           By default, lambda is 1 (full upwind stabilization)
!           ---------------------------------------------------
            lambdaStab = 1.0_RP

         end if
!
!        If central fluxes are used, set lambdaStab to zero
!        --------------------------------------------------
         if (whichRiemannSolver .eq. RIEMANN_CENTRAL) lambdaStab = 0.0_RP
!
!        ----------------------------
!        Set up an averaging function
!        ----------------------------
         if (controlVariables % containsKey(AVG_NAME_KEY)) then

            keyword = controlVariables % stringValueForKey(AVG_NAME_KEY, KEYWORD_LENGTH)
            call toLower(keyword)

            select case (keyword)
            case (STANDARD_AVG_NAME)
               AveragedStates => StandardAverage
               TwoPointFlux => StandardDG_TwoPointFlux
               whichAverage = STANDARD_AVG

            case (SKEWSYMMETRIC1_AVG_NAME)
               AveragedStates => SkewSymmetric1Average
               TwoPointFlux => SkewSymmetric1DG_TwoPointFlux
               whichAverage = SKEWSYMMETRIC1_AVG

            case (SKEWSYMMETRIC2_AVG_NAME)
               AveragedStates => SkewSymmetric2Average
               TwoPointFlux => SkewSymmetric2DG_TwoPointFlux
               whichAverage = SKEWSYMMETRIC2_AVG

            case default
               print*, "Averaging not recognized."
               errorMessage(STD_OUT)
               error stop
            end select

         else
!
!           Select standard by default
!           --------------------------
            AveragedStates => StandardAverage
            TwoPointFlux => StandardDG_TwoPointFlux
            whichAverage = STANDARD_AVG

         end if

      END SUBROUTINE SetRiemannSolver

      subroutine DescribeRiemannSolver
!
!        -------
!        Modules
!        -------
         use RiemannSolvers_iNSKeywordsModule


         select case (whichAverage)
         case (STANDARD_AVG)
            write(STD_OUT,'(30X,A,A30,A)') "->","Averaging function: ","Standard"

         case (SKEWSYMMETRIC1_AVG)
            write(STD_OUT,'(30X,A,A30,A)') "->","Averaging function: ","Skew-symmetric 1"

         case (SKEWSYMMETRIC2_AVG)
            write(STD_OUT,'(30X,A,A30,A)') "->","Averaging function: ","Skew-symmetric 2"

         end select

         select case (whichRiemannSolver)
         case (RIEMANN_CENTRAL)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Central"

         case (RIEMANN_LXF)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Lax-Friedrichs"

         case (RIEMANN_EXACT)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Exact"

         end select

         write(STD_OUT,'(30X,A,A30,F10.3)') "->","Lambda stabilization: ", lambdaStab

      end subroutine DescribeRiemannSolver
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
!        Riemann solvers
!        ---------------
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine CentralRiemannSolver(QLeft, QRight, nHat, t1, t2, flux)
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: f(1:NCONS), g(1:NCONS), h(1:NCONS)
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
!         rhoL = QLeft(INSRHO)
!         invRhoL = 1.0_RP / rhoL
!         uL = invRhoL * (QLeft(INSRHOU) * nHat(1) + QLeft(INSRHOV) * nHat(2) + QLeft(INSRHOW) * nHat(3))
!         vL = invRhoL * (QLeft(INSRHOU) * t1(1)   + QLeft(INSRHOV) * t1(2)   + QLeft(INSRHOW) * t1(3))
!         wL = invRhoL * (QLeft(INSRHOU) * t2(1)   + QLeft(INSRHOV) * t2(2)   + QLeft(INSRHOW) * t2(3))
!         pL = QLeft(INSP)
!
!         rhoR = QRight(INSRHO)
!         invRhoR = 1.0_RP / rhoR
!         uR = invRhoR * (QRight(INSRHOU) * nHat(1) + QRight(INSRHOV) * nHat(2) + QRight(INSRHOW) * nHat(3))
!         vR = invRhoR * (QRight(INSRHOU) * t1(1)   + QRight(INSRHOV) * t1(2)   + QRight(INSRHOW) * t1(3))
!         wR = invRhoR * (QRight(INSRHOU) * t2(1)   + QRight(INSRHOV) * t2(2)   + QRight(INSRHOW) * t2(3))
!         pR = QRight(INSP)
!!
!!        Perform the average using the averaging function
!!        ------------------------------------------------
!         QLRot = (/ rhoL, uL, vL, wL, pL /)
!         QRRot = (/ rhoR, uR, vR, wR, pR /)
!         call AveragedStates(QLRot, QRRot, pL, pR, rhoL, rhoR, flux)
!!
!!        ************************************************
!!        Return momentum equations to the cartesian frame
!!        ************************************************
!!
!         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

          call AveragedStates(QLeft, QRight, f, g, h)

          flux = f * nHat(1) + g*nHat(2) + h*nHat(3)

      end subroutine CentralRiemannSolver

      subroutine LxFRiemannSolver(QLeft, QRight, nHat, t1, t2, flux)
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: rhoL, uL, vL, wL, pL, invRhoL
         real(kind=RP)  :: rhoR, uR, vR, wR, pR, invRhoR
         real(kind=RP)  :: QLRot(NCONS), QRRot(NCONS)
         real(kind=RP)  :: stab(NCONS), lambdaMax
!!
!!        Rotate the variables to the face local frame using normal and tangent vectors
!!        -----------------------------------------------------------------------------
!         rhoL = QLeft(INSRHO)
!         invRhoL = 1.0_RP / rhoL
!         uL = invRhoL * (QLeft(INSRHOU) * nHat(1) + QLeft(INSRHOV) * nHat(2) + QLeft(INSRHOW) * nHat(3))
!         vL = invRhoL * (QLeft(INSRHOU) * t1(1)   + QLeft(INSRHOV) * t1(2)   + QLeft(INSRHOW) * t1(3))
!         wL = invRhoL * (QLeft(INSRHOU) * t2(1)   + QLeft(INSRHOV) * t2(2)   + QLeft(INSRHOW) * t2(3))
!         pL = QLeft(INSP)
!
!         rhoR = QRight(INSRHO)
!         invRhoR = 1.0_RP / rhoR
!         uR = invRhoR * (QRight(INSRHOU) * nHat(1) + QRight(INSRHOV) * nHat(2) + QRight(INSRHOW) * nHat(3))
!         vR = invRhoR * (QRight(INSRHOU) * t1(1)   + QRight(INSRHOV) * t1(2)   + QRight(INSRHOW) * t1(3))
!         wR = invRhoR * (QRight(INSRHOU) * t2(1)   + QRight(INSRHOV) * t2(2)   + QRight(INSRHOW) * t2(3))
!         pR = QRight(INSP)
!!
!!        Perform the average using the averaging function
!!        ------------------------------------------------
!         QLRot = (/ rhoL, uL, vL, wL, pL /)
!         QRRot = (/ rhoR, uR, vR, wR, pR /)
!         call AveragedStates(QLRot, QRRot, pL, pR, rhoL, rhoR, flux)
!!
!!        Compute the Lax-Friedrichs stabilization
!!        ----------------------------------------
!         lambdaMax = max(uL + sqrt(uL**2+4.0_RP*thermodynamics % rho0c02/rhoL), &
!                      uR + sqrt(uR**2+4.0_RP*thermodynamics % rho0c02/rhoR)    ) 
!
!         stab = 0.5_RP * lambdaMax * (QRRot - QLRot)
!
!         flux = flux - stab
!!
!!        ************************************************
!!        Return momentum equations to the cartesian frame
!!        ************************************************
!!
!         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

print*, "LxF Riemann solver not implemented"
error stop

      end subroutine LxFRiemannSolver

      subroutine ExactRiemannSolver(QLeft, QRight, nHat, t1, t2, flux)
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: rhoL, uL, vL, wL, pL, invRhoL, lambdaMinusL, lambdaPlusL
         real(kind=RP)  :: rhoR, uR, vR, wR, pR, invRhoR, lambdaMinusR, lambdaPlusR
         real(kind=RP)  :: rhoStarL, rhoStarR, uStar, pStar, rhoStar, vStar, wStar
         real(kind=RP)  :: QLRot(NCONS), QRRot(NCONS)
         real(kind=RP)  :: stab(NCONS), lambdaMax
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         rhoL = QLeft(INSRHO)
         invRhoL = 1.0_RP / rhoL
         uL = invRhoL * (QLeft(INSRHOU) * nHat(1) + QLeft(INSRHOV) * nHat(2) + QLeft(INSRHOW) * nHat(3))
         vL = invRhoL * (QLeft(INSRHOU) * t1(1)   + QLeft(INSRHOV) * t1(2)   + QLeft(INSRHOW) * t1(3))
         wL = invRhoL * (QLeft(INSRHOU) * t2(1)   + QLeft(INSRHOV) * t2(2)   + QLeft(INSRHOW) * t2(3))
         pL = QLeft(INSP)

         rhoR = QRight(INSRHO)
         invRhoR = 1.0_RP / rhoR
         uR = invRhoR * (QRight(INSRHOU) * nHat(1) + QRight(INSRHOV) * nHat(2) + QRight(INSRHOW) * nHat(3))
         vR = invRhoR * (QRight(INSRHOU) * t1(1)   + QRight(INSRHOV) * t1(2)   + QRight(INSRHOW) * t1(3))
         wR = invRhoR * (QRight(INSRHOU) * t2(1)   + QRight(INSRHOV) * t2(2)   + QRight(INSRHOW) * t2(3))
         pR = QRight(INSP)
!
!        Compute the Star Region
!        -----------------------
         lambdaMinusR = 0.5_RP * (uR - sqrt(uR*uR + 4.0_RP*thermodynamics % rho0c02/rhoR))
         lambdaPlusR  = 0.5_RP * (uR + sqrt(uR*uR + 4.0_RP*thermodynamics % rho0c02/rhoR))

         lambdaMinusL = 0.5_RP * (uL - sqrt(uL*uL + 4.0_RP*thermodynamics % rho0c02/rhoL))
         lambdaPlusL  = 0.5_RP * (uL + sqrt(uL*uL + 4.0_RP*thermodynamics % rho0c02/rhoL))

         uStar = (pR-pL+rhoR*uR*lambdaMinusR-rhoL*uL*lambdaPlusL)/(rhoR*lambdaMinusR - rhoL*lambdaPlusL)
         pStar = pR + rhoR*lambdaMinusR*(uR-uStar)
         rhoStarL = (rhoL*lambdaPlusL)/(uStar-lambdaMinusL)
         rhoStarR = (rhoR*lambdaMinusR)/(uStar - lambdaPlusR)

         if ( uStar .ge. 0.0_RP ) then
            rhoStar = rhoStarL
            vStar   = vL
            wStar   = wL

         else
            rhoStar = rhoStarR
            vStar   = vR
            wStar   = wR

         end if

         flux = [rhoStar*uStar, rhoStar*uStar*uStar + pStar, rhoStar*uStar*vStar, rhoStar*uStar*wStar, thermodynamics % rho0c02 * uStar]
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

      end subroutine ExactRiemannSolver
!
!////////////////////////////////////////////////////////////////////////////////////////////
!
!        Averaged states functions
!        -------------------------
!
!     To this averages, the states QLeft and QRight velocities have already been rotated,
!  so that only the first flux component has to be computed (Rotational invariance, see
!  E.F. Toro - Riemann solvers and numerical methods for fluid dynamics. Page 105).
!
!  Implemented two-point averages are:
!     -> Standard: Central fluxes
!     -> Skew-symmetric 1
!     -> Skew-symmetric 2
!
!////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine StandardAverage(QLeft, QRight, f, g, h)
!
!        *********************************************************************
!           Computes the standard average of the two states:
!              F* = {{F}} = 0.5 * (FL + FR)
!
!           State vectors are rotated.
!        *********************************************************************
!
         use Physics_iNS, only: iEulerXFlux
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(out)      :: f(1:NCONS), g(1:NCONS), h(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: fL(NCONS), fR(NCONS)
!
!        Compute the flux
!        ----------------
!         fL(INSRHO)  = QLeft(INSRHO) * QLeft(INSRHOU)
!         fL(INSRHOU) = fL(INSRHO) * QLeft(INSRHOU) + QLeft(INSP)
!         fL(INSRHOV) = fL(INSRHO) * QLeft(INSRHOV)
!         fL(INSRHOW) = fL(INSRHO) * QLeft(INSRHOW)
!         fL(INSP)    = thermodynamics % rho0c02 * QLeft(INSP)
!
!         fR(INSRHO)  = QRight(INSRHO) * QRight(INSRHOU)
!         fR(INSRHOU) = fR(INSRHO) * QRight(INSRHOU) + QRight(INSP)
!         fR(INSRHOV) = fR(INSRHO) * QRight(INSRHOV)
!         fR(INSRHOW) = fR(INSRHO) * QRight(INSRHOW)
!         fR(INSP)    = thermodynamics % rho0c02 * QRight(INSP)
!
!         flux = 0.5_RP * (fL + fR)

         print*, "Standard average not implemented"
         error stop

      end subroutine StandardAverage

      subroutine SkewSymmetric1Average(QLeft, QRight, f, g, h)
!
!        *********************************************************************
!        *********************************************************************
!
         use Physics_iNS, only: iEulerXFlux
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(out)      :: f(1:NCONS), g(1:NCONS), h(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: rhou, u

!         rhou = 0.5_RP * (QLeft(INSRHO) * QLeft(INSRHOU) + QRight(INSRHO) * QRight(INSRHOU))
!         u    = 0.5_RP * (QLeft(INSRHOU) + QRight(INSRHOU))
!!
!!        Compute the flux
!!        ----------------
!         flux(INSRHO) = rhou
!         flux(INSRHOU) = rhou * u + 0.5_RP * (QLeft(INSP) + QRight(INSP))
!         flux(INSRHOU) = rhou * 0.5_RP * (QLeft(INSRHOV) + QRight(INSRHOV))
!         flux(INSRHOU) = rhou * 0.5_RP * (QLeft(INSRHOW) + QRight(INSRHOW))
!         flux(INSP)    = thermodynamics % rho0c02 * u
!
         print*, "SkewSymmetric 1 average not implemented"
         error stop
      end subroutine SkewSymmetric1Average

      subroutine SkewSymmetric2Average(QLeft, QRight, f, g, h)
!
!        *********************************************************************
!        *********************************************************************
!
         use Physics_iNS, only: iEulerXFlux
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(out)      :: f(1:NCONS), g(1:NCONS), h(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: invRhoL, invRhoR
         real(kind=RP)  :: rho, u, v, w, p

         invRhoL = 1.0_RP / QLeft(INSRHO)    ; invRhoR = 1.0_RP / QRight(INSRHO)

         rho = 0.5_RP * (QLeft(INSRHO) + QRight(INSRHO))
         u   = 0.5_RP * (invRhoL * QLeft(INSRHOU) + invRhoR * QRight(INSRHOU))
         v   = 0.5_RP * (invRhoL * QLeft(INSRHOV) + invRhoR * QRight(INSRHOV))
         w   = 0.5_RP * (invRhoL * QLeft(INSRHOW) + invRhoR * QRight(INSRHOW))
         p   = 0.5_RP * (QLeft(INSP) + QRight(INSP))

         f(INSRHO) = rho*u
         f(INSRHOU) = f(INSRHO)*u+p
         f(INSRHOV) = f(INSRHO)*v
         f(INSRHOW) = f(INSRHO)*w
         f(INSP)    = thermodynamics % rho0c02 * u

         g(INSRHO) = rho*v
         g(INSRHOU) = g(INSRHO)*u
         g(INSRHOV) = g(INSRHO)*v+p
         g(INSRHOW) = g(INSRHO)*w
         g(INSP)    = thermodynamics % rho0c02 * v

         h(INSRHO) = rho*w
         h(INSRHOU) = h(INSRHO)*u
         h(INSRHOV) = h(INSRHO)*v
         h(INSRHOW) = h(INSRHO)*w + p
         h(INSP)    = thermodynamics % rho0c02 * w


         

!         u    = 0.5_RP * (QLeft(INSRHOU) + QRight(INSRHOU))
!!
!!        Compute the flux
!!        ----------------
!         flux(INSRHO) = 0.5_RP * (QLeft(INSRHO) + QRight(INSRHO)) * u
!         flux(INSRHOU) = flux(INSRHOU) * u + 0.5_RP * (QLeft(INSP) + QRight(INSP))
!         flux(INSRHOU) = flux(INSRHOU) * 0.5_RP * (QLeft(INSRHOV) + QRight(INSRHOV))
!         flux(INSRHOU) = flux(INSRHOU) * 0.5_RP * (QLeft(INSRHOW) + QRight(INSRHOW))
!         flux(INSP)    = thermodynamics % rho0c02 * u

      end subroutine SkewSymmetric2Average
!
!///////////////////////////////////////////////////////////////////////
!
!        Volumetric two-point fluxes
!        ---------------------------
!///////////////////////////////////////////////////////////////////////
!
      subroutine StandardDG_TwoPointFlux(QL,QR,JaL,JaR, fSharp)
         use SMConstants
         use PhysicsStorage_iNS
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: rhoL, uL, vL, wL, pL, invRhoL
         real(kind=RP)     :: rhoR, uR, vR, wR, pR, invRhoR
         real(kind=RP)     :: Ja(1:NDIM)
         real(kind=RP)     :: f(NCONS), g(NCONS), h(NCONS)

         rhoL    = QL(INSRHO)
         rhoR    = QR(INSRHO)
         invRhoL = 1.0_RP / rhoL         ; invRhoR = 1.0_RP / rhoR
         uL      = QL(INSRHOU) * invRhoL ; uR      = QR(INSRHOU) * invRhoR
         vL      = QL(INSRHOV) * invRhoL ; vR      = QR(INSRHOV) * invRhoR
         wL      = QL(INSRHOW) * invRhoL ; wR      = QR(INSRHOW) * invRhoR
         pL      = QL(INSP)              ; pR      = QR(INSP)

!
!        Average metrics: (Note: Here all average (1/2)s are accounted later)
!        ---------------
         Ja = (JaL + JaR)
!
!        Compute the flux
!        ----------------
         f(INSRHO)  = rhoL*uL         + rhoR*uR
         f(INSRHOU) = rhoL*uL*uL + pL + rhoR*uR*uR + pR
         f(INSRHOV) = rhoL*uL*vL      + rhoR*uR*vR
         f(INSRHOW) = rhoL*uL*wL      + rhoR*uR*wR
         f(INSP)    = thermodynamics % rho0c02 * (uL + uR)

         g(INSRHO)  = rhoL*vL         + rhoR*vR
         g(INSRHOU) = rhoL*vL*uL      + rhoR*vR*uR
         g(INSRHOV) = rhoL*vL*vL + pL + rhoR*vR*vR + pR
         g(INSRHOW) = rhoL*vL*wL      + rhoR*vR*wR
         g(INSP)    = thermodynamics % rho0c02 * (vL + vR)

         h(INSRHO)  = rhoL*wL         + rhoR*wR
         h(INSRHOU) = rhoL*wL*uL      + rhoR*wR*uR
         h(INSRHOV) = rhoL*wL*vL      + rhoR*wR*vR
         h(INSRHOW) = rhoL*wL*wL + pL + rhoR*wR*wR + pR
         h(INSP)    = thermodynamics % rho0c02 * (wL + wR)
!
!        Compute the sharp flux: (And account for the (1/2)^2)
!        ----------------------
         fSharp = 0.25_RP * ( f*Ja(IX) + g*Ja(IY) + h*Ja(IZ) )

      end subroutine StandardDG_TwoPointFlux

      subroutine SkewSymmetric1DG_TwoPointFlux(QL,QR,JaL,JaR, fSharp)
         use SMConstants
         use PhysicsStorage_iNS
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: rhoL, uL, vL, wL, pL, invRhoL
         real(kind=RP) :: rhoR, uR, vR, wR, pR, invRhoR
         real(kind=RP) :: rho, u, v, w, p, rhou, rhov, rhow
         real(kind=RP) :: Ja(1:NDIM)
         real(kind=RP) :: f(NCONS), g(NCONS), h(NCONS)

         rhoL    = QL(INSRHO)
         rhoR    = QR(INSRHO)
         invRhoL = 1.0_RP / rhoL         ; invRhoR = 1.0_RP / rhoR
         uL      = QL(INSRHOU) * invRhoL ; uR      = QR(INSRHOU) * invRhoR
         vL      = QL(INSRHOV) * invRhoL ; vR      = QR(INSRHOV) * invRhoR
         wL      = QL(INSRHOW) * invRhoL ; wR      = QR(INSRHOW) * invRhoR
         pL      = QL(INSP)              ; pR      = QR(INSP)

         rho = 0.5_RP * (rhoL + rhoR)
         u   = 0.5_RP * (uL + uR)
         v   = 0.5_RP * (vL + vR)
         w   = 0.5_RP * (wL + wR)
         p   = 0.5_RP * (pL + pR)

         rhou = 0.5_RP * (QL(INSRHOU) + QR(INSRHOU))
         rhov = 0.5_RP * (QL(INSRHOV) + QR(INSRHOV))
         rhow = 0.5_RP * (QL(INSRHOW) + QR(INSRHOW))
!
!        Average metrics
!        ---------------
         Ja = 0.5_RP * (JaL + JaR)
!
!        Compute the flux
!        ----------------
         f(INSRHO)  = rhou
         f(INSRHOU) = rhou*u + p
         f(INSRHOV) = rhou*v
         f(INSRHOW) = rhou*w
         f(INSP)    = thermodynamics % rho0c02 * u

         g(INSRHO)  = rhov
         g(INSRHOU) = rhov*u
         g(INSRHOV) = rhov*v + p
         g(INSRHOW) = rhov*w
         g(INSP)    = thermodynamics % rho0c02 * v

         h(INSRHO)  = rhow
         h(INSRHOU) = rhow*u
         h(INSRHOV) = rhow*v
         h(INSRHOW) = rhow*w + p
         h(INSP)    = thermodynamics % rho0c02 * w
!
!        Compute the sharp flux: (And account for the (1/2)^2)
!        ----------------------
         fSharp = f*Ja(IX) + g*Ja(IY) + h*Ja(IZ)

      end subroutine SkewSymmetric1DG_TwoPointFlux

      subroutine SkewSymmetric2DG_TwoPointFlux(QL,QR,JaL,JaR, fSharp)
         use SMConstants
         use PhysicsStorage_iNS
         implicit none
         real(kind=RP), intent(in)       :: QL(1:NCONS)
         real(kind=RP), intent(in)       :: QR(1:NCONS)
         real(kind=RP), intent(in)       :: JaL(1:NDIM)
         real(kind=RP), intent(in)       :: JaR(1:NDIM)
         real(kind=RP), intent(out)      :: fSharp(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: rhoL, uL, vL, wL, pL, invRhoL
         real(kind=RP) :: rhoR, uR, vR, wR, pR, invRhoR
         real(kind=RP) :: rho, u, v, w, p
         real(kind=RP) :: Ja(1:NDIM)
         real(kind=RP) :: f(NCONS), g(NCONS), h(NCONS)

         rhoL    = QL(INSRHO)
         rhoR    = QR(INSRHO)

         invRhoL = 1.0_RP / rhoL         ; invRhoR = 1.0_RP / rhoR
         uL      = QL(INSRHOU) * invRhoL ; uR      = QR(INSRHOU) * invRhoR
         vL      = QL(INSRHOV) * invRhoL ; vR      = QR(INSRHOV) * invRhoR
         wL      = QL(INSRHOW) * invRhoL ; wR      = QR(INSRHOW) * invRhoR
         pL      = QL(INSP)              ; pR      = QR(INSP)

         rho = 0.5_RP * (rhoL + rhoR)
         u   = 0.5_RP * (uL + uR)
         v   = 0.5_RP * (vL + vR)
         w   = 0.5_RP * (wL + wR)
         p   = 0.5_RP * (pL + pR)
!
!        Average metrics
!        ---------------
         Ja = 0.5_RP * (JaL + JaR)
!
!        Compute the flux
!        ----------------
         f(INSRHO)  = rho*u
         f(INSRHOU) = rho*u*u + p
         f(INSRHOV) = rho*u*v
         f(INSRHOW) = rho*u*w
         f(INSP)    = thermodynamics % rho0c02 * u

         g(INSRHO)  = rho*v
         g(INSRHOU) = rho*v*u
         g(INSRHOV) = rho*v*v + p
         g(INSRHOW) = rho*v*w
         g(INSP)    = thermodynamics % rho0c02 * v

         h(INSRHO)  = rho*w
         h(INSRHOU) = rho*w*u
         h(INSRHOV) = rho*w*v
         h(INSRHOW) = rho*w*w + p
         h(INSP)    = thermodynamics % rho0c02 * w
!
!        Compute the sharp flux: (And account for the (1/2)^2)
!        ----------------------
         fSharp = f*Ja(IX) + g*Ja(IY) + h*Ja(IZ)

      end subroutine SkewSymmetric2DG_TwoPointFlux

end module RiemannSolvers_iNS