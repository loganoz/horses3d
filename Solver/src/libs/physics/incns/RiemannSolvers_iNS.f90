!
!//////////////////////////////////////////////////////
!
!   @File:    RiemannSolvers_iNS.f90
!   @Author:  Juan Manzanero (j.manzanero1992@gmail.com)
!   @Created: Tue Jun 19 17:39:26 2018
!   @Last revision date: Fri Jun 22 12:58:37 2018
!   @Last revision author: Juan Manzanero (j.manzanero1992@gmail.com)
!   @Last revision commit: 5fcd2e67947be854342011924f2897ed668cf53a
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
         real(kind=RP)  :: rhoL, uL, vL, wL, pL
         real(kind=RP)  :: rhoR, uR, vR, wR, pR
         real(kind=RP)  :: QLRot(NINC), QRRot(NINC)
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         rhoL = QLeft(INSRHO)
         uL = QLeft(INSU) * nHat(1) + QLeft(INSV) * nHat(2) + QLeft(INSW) * nHat(3)
         vL = QLeft(INSU) * t1(1)   + QLeft(INSV) * t1(2)   + QLeft(INSW) * t1(3)
         wL = QLeft(INSU) * t2(1)   + QLeft(INSV) * t2(2)   + QLeft(INSW) * t2(3)
         pL = QLeft(INSP)

         rhoR = QRight(INSRHO)
         uR = QRight(INSU) * nHat(1) + QRight(INSV) * nHat(2) + QRight(INSW) * nHat(3)
         vR = QRight(INSU) * t1(1)   + QRight(INSV) * t1(2)   + QRight(INSW) * t1(3)
         wR = QRight(INSU) * t2(1)   + QRight(INSV) * t2(2)   + QRight(INSW) * t2(3)
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
         real(kind=RP)  :: rhoL, uL, vL, wL, pL
         real(kind=RP)  :: rhoR, uR, vR, wR, pR
         real(kind=RP)  :: QLRot(NINC), QRRot(NINC)
         real(kind=RP)  :: stab(NINC), lambdaMax
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         rhoL = QLeft(INSRHO)
         uL = QLeft(INSU) * nHat(1) + QLeft(INSV) * nHat(2) + QLeft(INSW) * nHat(3)
         vL = QLeft(INSU) * t1(1)   + QLeft(INSV) * t1(2)   + QLeft(INSW) * t1(3)
         wL = QLeft(INSU) * t2(1)   + QLeft(INSV) * t2(2)   + QLeft(INSW) * t2(3)
         pL = QLeft(INSP)

         rhoR = QRight(INSRHO)
         uR = QRight(INSU) * nHat(1) + QRight(INSV) * nHat(2) + QRight(INSW) * nHat(3)
         vR = QRight(INSU) * t1(1)   + QRight(INSV) * t1(2)   + QRight(INSW) * t1(3)
         wR = QRight(INSU) * t2(1)   + QRight(INSV) * t2(2)   + QRight(INSW) * t2(3)
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
         real(kind=RP)  :: rhoL, uL, vL, wL, pL, lambdaMinusL, lambdaPlusL
         real(kind=RP)  :: rhoR, uR, vR, wR, pR, lambdaMinusR, lambdaPlusR
         real(kind=RP)  :: rhoStarL, rhoStarR, uStar, pStar
         real(kind=RP)  :: QLRot(NINC), QRRot(NINC), QStar(NINC)
         real(kind=RP)  :: stab(NINC), lambdaMax
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         rhoL = QLeft(INSRHO)
         uL = QLeft(INSU) * nHat(1) + QLeft(INSV) * nHat(2) + QLeft(INSW) * nHat(3)
         vL = QLeft(INSU) * t1(1)   + QLeft(INSV) * t1(2)   + QLeft(INSW) * t1(3)
         wL = QLeft(INSU) * t2(1)   + QLeft(INSV) * t2(2)   + QLeft(INSW) * t2(3)
         pL = QLeft(INSP)

         rhoR = QRight(INSRHO)
         uR = QRight(INSU) * nHat(1) + QRight(INSV) * nHat(2) + QRight(INSW) * nHat(3)
         vR = QRight(INSU) * t1(1)   + QRight(INSV) * t1(2)   + QRight(INSW) * t1(3)
         wR = QRight(INSU) * t2(1)   + QRight(INSV) * t2(2)   + QRight(INSW) * t2(3)
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
            Qstar = [rhoStarL, uStar, vL, wL, pStar]

         else
            QStar = [rhoStarR, uStar, vR, wR, pStar]

         end if

         flux = [QStar(INSRHO)*QStar(INSU),QStar(INSRHO)*QStar(INSU)*QStar(INSU) + QStar(INSP), &
                 QStar(INSRHO)*QStar(INSU)*QStar(INSV), QStar(INSRHO)*QStar(INSU)*QStar(INSW),  &
                 thermodynamics % rho0c02*QStar(INSU)] 
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
         flux(INSRHO) = 0.5_RP * (QLeft(INSRHO)*QLeft(INSU)                     + QRight(INSRHO)*QRight(INSU))
         flux(INSU)   = 0.5_RP * (QLeft(INSRHO)*POW2(QLeft(INSU)) + QLeft(INSP) + QRight(INSRHO)*POW2(QRight(INSU)) + QRight(INSP))
         flux(INSV)   = 0.5_RP * (QLeft(INSRHO)*QLeft(INSU)*QLeft(INSV)         + QRight(INSRHO)*QRight(INSU)*QRight(INSV))
         flux(INSW)   = 0.5_RP * (QLeft(INSRHO)*QLeft(INSU)*QLeft(INSW)         + QRight(INSRHO)*QRight(INSU)*QRight(INSW))
         flux(INSP)   = 0.5_RP * (QLeft(INSU)                                   + QRight(INSU)) * thermodynamics % rho0c02

      end subroutine StandardAverage

end module RiemannSolvers_iNS
