!
!//////////////////////////////////////////////////////
!
!   @File:    RiemannSolvers_iNS.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Jun 19 17:39:26 2018
!   @Last revision date: Wed Jul 18 10:33:20 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 4977ebc1252872ccf3ec1e535ebb8619da12e2c8
!
!//////////////////////////////////////////////////////
!
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
         real(kind=RP), intent(in)       :: QLeft(1:NINC)
         real(kind=RP), intent(in)       :: QRight(1:NINC)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NINC)
      end subroutine RiemannSolverFCN

      subroutine AveragedStatesFCN(QLeft, QRight, pL, pR, invRhoL, invRhoR, flux) 
         use SMConstants
         use PhysicsStorage_iNS
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NINC)
         real(kind=RP), intent(in)       :: QRight(1:NINC)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NINC)
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

         case (SKEWSYMMETRIC_SPLIT)
            AveragedStates => StandardAverage
            whichAverage = STANDARD_SPLIT


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
         real(kind=RP), intent(in)       :: QLeft(1:NINC)
         real(kind=RP), intent(in)       :: QRight(1:NINC)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NINC)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: rhoL, uL, vL, wL, pL, invRhoL
         real(kind=RP)  :: rhoR, uR, vR, wR, pR, invRhoR
         real(kind=RP)  :: QLRot(NINC), QRRot(NINC)
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
!        Perform the average using the averaging function
!        ------------------------------------------------
         QLRot = (/ rhoL, uL, vL, wL, pL /)
         QRRot = (/ rhoR, uR, vR, wR, pR /)
         call AveragedStates(QLRot, QRRot, pL, pR, rhoL, rhoR, flux)
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

      end subroutine CentralRiemannSolver

      subroutine LxFRiemannSolver(QLeft, QRight, nHat, t1, t2, flux)
         implicit none 
         real(kind=RP), intent(in)       :: QLeft(1:NINC)
         real(kind=RP), intent(in)       :: QRight(1:NINC)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NINC)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: rhoL, uL, vL, wL, pL, invRhoL
         real(kind=RP)  :: rhoR, uR, vR, wR, pR, invRhoR
         real(kind=RP)  :: QLRot(NINC), QRRot(NINC)
         real(kind=RP)  :: stab(NINC), lambdaMax
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
!        Perform the average using the averaging function
!        ------------------------------------------------
         QLRot = (/ rhoL, uL, vL, wL, pL /)
         QRRot = (/ rhoR, uR, vR, wR, pR /)
         call AveragedStates(QLRot, QRRot, pL, pR, rhoL, rhoR, flux)
!
!        Compute the Lax-Friedrichs stabilization
!        ----------------------------------------
         lambdaMax = max(uL + sqrt(uL**2+4.0_RP*thermodynamics % rho0c02/rhoL), &
                      uR + sqrt(uR**2+4.0_RP*thermodynamics % rho0c02/rhoR)    ) 

         stab = 0.5_RP * lambdaMax * (QRRot - QLRot)

         flux = flux - stab
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

      end subroutine LxFRiemannSolver

      subroutine ExactRiemannSolver(QLeft, QRight, nHat, t1, t2, flux)
         implicit none 
         real(kind=RP), intent(in)       :: QLeft(1:NINC)
         real(kind=RP), intent(in)       :: QRight(1:NINC)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NINC)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: rhoL, uL, vL, wL, pL, invRhoL, lambdaMinusL, lambdaPlusL
         real(kind=RP)  :: rhoR, uR, vR, wR, pR, invRhoR, lambdaMinusR, lambdaPlusR
         real(kind=RP)  :: rhoStarL, rhoStarR, uStar, pStar, rhoStar, vStar, wStar
         real(kind=RP)  :: QLRot(NINC), QRRot(NINC) 
         real(kind=RP)  :: stab(NINC), lambdaMax
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
      subroutine StandardAverage(QLeft, QRight, pL, pR, invRhoL, invRhoR, flux) 
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
         real(kind=RP), intent(in)       :: QLeft(1:NINC)
         real(kind=RP), intent(in)       :: QRight(1:NINC)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NINC)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: fL(NINC), fR(NINC)
!
!        Compute the flux
!        ----------------
         call iEulerXFlux(QLeft , fL)
         call iEulerXFlux(QRight, fR)
      
         flux = 0.5_RP * (fL + fR)

      end subroutine StandardAverage

end module RiemannSolvers_iNS
