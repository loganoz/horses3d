!
!//////////////////////////////////////////////////////
!
!   @File:    RiemannSolvers_MU.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Jun 19 17:39:26 2018
!   @Last revision date: Wed Jul 18 10:33:20 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 4977ebc1252872ccf3ec1e535ebb8619da12e2c8
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module RiemannSolvers_MU
   use SMConstants
   use Physics_MU
   use PhysicsStorage_MU
   use VariableConversion_MU
   use FluidData_MU

   private 
   public RiemannSolver, SetRiemannSolver, ExactRiemannSolver

   abstract interface
      subroutine RiemannSolverFCN(QLeft, QRight, rhoL, rhoR, muL, muR, nHat, t1, t2, fL,fR)
         use SMConstants
         use PhysicsStorage_MU
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: rhoL
         real(kind=RP), intent(in)       :: rhoR
         real(kind=RP), intent(in)       :: muL
         real(kind=RP), intent(in)       :: muR
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: fL(1:NCONS)
         real(kind=RP), intent(out)      :: fR(1:NCONS)
      end subroutine RiemannSolverFCN
   end interface

   procedure(RiemannSolverFCN)     , pointer  :: RiemannSolver      => NULL()

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

         case(RIEMANN_EXACT)
            RiemannSolver => ExactRiemannSolver

         case default
            print*, "Undefined choice of Riemann Solver."
            print*, "Options available are:"
            print*, "   * Central"
            print*, "   * Exact"
            errorMessage(STD_OUT)
            STOP
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
      subroutine CentralRiemannSolver(QLeft, QRight, rhoL, rhoR, muL, muR, nHat, t1, t2, fL, fR)
         implicit none 
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: rhoL
         real(kind=RP), intent(in)       :: rhoR
         real(kind=RP), intent(in)       :: muL
         real(kind=RP), intent(in)       :: muR
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: fL(1:NCONS)
         real(kind=RP), intent(out)      :: fR(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: cL, uL, vL, wL, pL, invSqrtRhoL
         real(kind=RP) :: cR, uR, vR, wR, pR, invSqrtRhoR
!
!        Left state variables and fluxes
!        -------------------------------
         invSqrtRhoL = 1.0_RP / sqrt(rhoL)
         cL = QLeft(IMC)
         uL = invSqrtRhoL * (QLeft(IMSQRHOU) * nHat(1) + QLeft(IMSQRHOV) * nHat(2) + QLeft(IMSQRHOW) * nHat(3))
         vL = invSqrtRhoL * (QLeft(IMSQRHOU) * t1(1)   + QLeft(IMSQRHOV) * t1(2)   + QLeft(IMSQRHOW) * t1(3))
         wL = invSqrtRhoL * (QLeft(IMSQRHOU) * t2(1)   + QLeft(IMSQRHOV) * t2(2)   + QLeft(IMSQRHOW) * t2(3))
         pL = QLeft(IMP)

         fL(IMC)      = uL*cL
         fL(IMSQRHOU) = 0.5_RP*rhoL*uL*uL + pL
         fL(IMSQRHOV) = 0.5_RP*rhoL*uL*vL
         fL(IMSQRHOW) = 0.5_RP*rhoL*uL*wL
         fL(IMP)      = thermodynamics % rho0c02 * uL

!
!        Right state variables and fluxes
!        --------------------------------
         invSqrtRhoR = 1.0_RP / sqrt(rhoR)
         cR = QRight(IMC)
         uR = invSqrtRhoR * (QRight(IMSQRHOU) * nHat(1) + QRight(IMSQRHOV) * nHat(2) + QRight(IMSQRHOW) * nHat(3))
         vR = invSqrtRhoR * (QRight(IMSQRHOU) * t1(1)   + QRight(IMSQRHOV) * t1(2)   + QRight(IMSQRHOW) * t1(3))
         wR = invSqrtRhoR * (QRight(IMSQRHOU) * t2(1)   + QRight(IMSQRHOV) * t2(2)   + QRight(IMSQRHOW) * t2(3))
         pR = QRight(IMP)

         fR(IMC)      = uR*cR
         fR(IMSQRHOU) = 0.5_RP*rhoR*uR*uR + pR
         fR(IMSQRHOV) = 0.5_RP*rhoR*uR*vR
         fR(IMSQRHOW) = 0.5_RP*rhoR*uR*wR
         fR(IMP)      = thermodynamics % rho0c02 * uR
!
!        Perform the average and rotation
!        --------------------------------
         fL = 0.5_RP*(fL + fR)

         fR = fL
!
!        Add the non-conservative term
!        -----------------------------          
         fL(IMSQRHOU) = fL(IMSQRHOU) + 0.5_RP*cL*(muR-muL) + 0.25_RP*rhoL*uL*(uR-uL)
         fL(IMSQRHOV) = fL(IMSQRHOV) + 0.25_RP*rhoL*uL*(vR-vL)
         fL(IMSQRHOW) = fL(IMSQRHOW) + 0.25_RP*rhoL*uL*(wR-wL)

         fR(IMSQRHOU) = fR(IMSQRHOU) + 0.5_RP*cR*(muL-muR) + 0.25_RP*rhoR*uR*(uL-uR)
         fR(IMSQRHOV) = fR(IMSQRHOV) + 0.25_RP*rhoR*uR*(vL-vR)
         fR(IMSQRHOW) = fR(IMSQRHOW) + 0.25_RP*rhoR*uR*(wL-wR)

         fL(IMSQRHOU:IMSQRHOW) = nHat*fL(IMSQRHOU) + t1*fL(IMSQRHOV) + t2*fL(IMSQRHOW)
         fR(IMSQRHOU:IMSQRHOW) = nHat*fR(IMSQRHOU) + t1*fR(IMSQRHOV) + t2*fR(IMSQRHOW)

      end subroutine CentralRiemannSolver

      subroutine ExactRiemannSolver(QLeft, QRight, rhoL, rhoR, muL, muR, nHat, t1, t2, fL, fR)
         implicit none 
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: rhoL
         real(kind=RP), intent(in)       :: rhoR
         real(kind=RP), intent(in)       :: muL
         real(kind=RP), intent(in)       :: muR
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: fL(1:NCONS)
         real(kind=RP), intent(out)      :: fR(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: uL, vL, wL, pL, invRhoL, lambdaMinusL, lambdaPlusL
         real(kind=RP)  :: uR, vR, wR, pR, invRhoR, lambdaMinusR, lambdaPlusR
         real(kind=RP)  :: rhoStarL, rhoStarR, uStar, pStar, rhoStar, vStar, wStar
         real(kind=RP)  :: QLRot(NCONS), QRRot(NCONS) 
         real(kind=RP)  :: stab(NCONS), lambdaMax
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
!         invRhoL = 1.0_RP / rhoL
!         uL = invRhoL * (QLeft(INSRHOU) * nHat(1) + QLeft(INSRHOV) * nHat(2) + QLeft(INSRHOW) * nHat(3))
!         vL = invRhoL * (QLeft(INSRHOU) * t1(1)   + QLeft(INSRHOV) * t1(2)   + QLeft(INSRHOW) * t1(3))
!         wL = invRhoL * (QLeft(INSRHOU) * t2(1)   + QLeft(INSRHOV) * t2(2)   + QLeft(INSRHOW) * t2(3))
!         pL = QLeft(INSP)
!
!         invRhoR = 1.0_RP / rhoR
!         uR = invRhoR * (QRight(INSRHOU) * nHat(1) + QRight(INSRHOV) * nHat(2) + QRight(INSRHOW) * nHat(3))
!         vR = invRhoR * (QRight(INSRHOU) * t1(1)   + QRight(INSRHOV) * t1(2)   + QRight(INSRHOW) * t1(3))
!         wR = invRhoR * (QRight(INSRHOU) * t2(1)   + QRight(INSRHOV) * t2(2)   + QRight(INSRHOW) * t2(3))
!         pR = QRight(INSP)
!!
!!        Compute the Star Region
!!        -----------------------
!         lambdaMinusR = 0.5_RP * (uR - sqrt(uR*uR + 4.0_RP*thermodynamics % rho0c02/rhoR))
!         lambdaPlusR  = 0.5_RP * (uR + sqrt(uR*uR + 4.0_RP*thermodynamics % rho0c02/rhoR))
!
!         lambdaMinusL = 0.5_RP * (uL - sqrt(uL*uL + 4.0_RP*thermodynamics % rho0c02/rhoL))
!         lambdaPlusL  = 0.5_RP * (uL + sqrt(uL*uL + 4.0_RP*thermodynamics % rho0c02/rhoL))
!
!         uStar = (pR-pL+rhoR*uR*lambdaMinusR-rhoL*uL*lambdaPlusL)/(rhoR*lambdaMinusR - rhoL*lambdaPlusL)
!         pStar = pR + rhoR*lambdaMinusR*(uR-uStar)
!         rhoStarL = (rhoL*lambdaPlusL)/(uStar-lambdaMinusL)
!         rhoStarR = (rhoR*lambdaMinusR)/(uStar - lambdaPlusR)
!
!         if ( uStar .ge. 0.0_RP ) then
!            rhoStar = rhoStarL
!            vStar   = vL
!            wStar   = wL
!
!         else
!            rhoStar = rhoStarR
!            vStar   = vR
!            wStar   = wR
!
!         end if
!
!         flux = [rhoStar*uStar, rhoStar*uStar*uStar + pStar, rhoStar*uStar*vStar, rhoStar*uStar*wStar, thermodynamics % rho0c02 * uStar]
!!
!!        ************************************************
!!        Return momentum equations to the cartesian frame
!!        ************************************************
!!
!         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

      end subroutine ExactRiemannSolver

end module RiemannSolvers_MU
