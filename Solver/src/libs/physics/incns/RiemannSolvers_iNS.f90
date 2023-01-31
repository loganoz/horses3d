#include "Includes.h"
module RiemannSolvers_iNS
   use SMConstants
   use Physics_iNS
   use PhysicsStorage_iNS
   use VariableConversion_iNS
   use FluidData_iNS

   private 
   public RiemannSolver, SetRiemannSolver, ExactRiemannSolver

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
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(out)      :: f(1:NCONS), g(1:NCONS), h(1:NCONS)
      end subroutine AveragedStatesFCN
   end interface

   procedure(RiemannSolverFCN)     , pointer  :: RiemannSolver      => NULL()
   procedure(AveragedStatesFCN)    , pointer  :: AveragedStates     => NULL()

   contains
      SUBROUTINE SetRiemannSolver(which, splitType)
!
!        **************************************************************
!              This subroutine is to set which Riemann solver is used.
!           the user cannot decide amongst the averaging function.
!           It is automatically selected depending on which split
!           form is enabled.
!           The user can choose the dissipation type:
!              None (central), Roe, Lax-Friedrichs, Rusanov
!
!           And the dissipation intensity, with the lambda stabilization
!           parameter. By default it is set to 1 (whole dissipation),
!           instead for central fluxes, which is 0 (no dissipation).
!        **************************************************************
!
         IMPLICIT NONE
         integer, intent(in) :: which
         integer, intent(in) :: splitType
         
         
         select case ( which )
         case(RIEMANN_CENTRAL)
            RiemannSolver => CentralRiemannSolver

         case(RIEMANN_LXF)
            RiemannSolver => LxFRiemannSolver

         case(RIEMANN_EXACT)
            RiemannSolver => ExactRiemannSolver

         case default
            print*, "Undefined choice of Riemann Solver."
            print*, "Options available are:"
            print*, "   * Central"
            print*, "   * Lax-Friedrichs"
            print*, "   * Exact"
            errorMessage(STD_OUT)
            STOP
         end select
!
!        Set up an averaging function
!        ----------------------------
         select case ( splitType )
         case (STANDARD_SPLIT)
            AveragedStates => StandardAverage
            whichAverage = STANDARD_SPLIT

         case (SKEWSYMMETRIC1_SPLIT)
            AveragedStates => SkewSymmetric1Average
            whichAverage = SKEWSYMMETRIC1_SPLIT

         case (SKEWSYMMETRIC2_SPLIT)
            AveragedStates => SkewSymmetric2Average
            whichAverage = SKEWSYMMETRIC2_SPLIT

         case default
            print*, "Split form not recognized"
            errorMessage(STD_OUT)
            stop
         end select
      
      END SUBROUTINE SetRiemannSolver
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
stop

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
!     -> Ducros
!     -> Morinishi
!     -> Kennedy and Gruber
!     -> Pirozzoli
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
         stop

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
         stop
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

end module RiemannSolvers_iNS