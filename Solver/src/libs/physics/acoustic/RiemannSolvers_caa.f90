#include "Includes.h"
module RiemannSolvers_CAAKeywordsModule

     integer,                       parameter :: KEYWORD_LENGTH           = 132
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_SOLVER_NAME_KEY  = "riemann solver"
     character(len=KEYWORD_LENGTH), parameter :: LAMBDA_STABILIZATION_KEY = "lambda stabilization"
!
!    --------------------------
!    Riemann solver definitions
!    --------------------------
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_CENTRAL_NAME    = "central"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_ROE_NAME        = "roe"
     ! character(len=KEYWORD_LENGTH), parameter :: RIEMANN_RUSANOV_NAME    = "rusanov"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_LXF_NAME        = "lax-friedrichs"

     enum, bind(C)
        ! enumerator :: RIEMANN_ROE = 1, RIEMANN_LXF, RIEMANN_RUSANOV
        enumerator :: RIEMANN_ROE = 1, RIEMANN_LXF
        enumerator :: RIEMANN_CENTRAL 
     end enum

end module RiemannSolvers_CAAKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
!
module RiemannSolvers_CAA
   use SMConstants
   use Physics_CAA
   use PhysicsStorage_CAA
   use VariableConversion_CAA
   use FluidData_CAA

   implicit none

   private
   public whichRiemannSolver
   public SetRiemannSolver, DescribeRiemannSolver
   public RiemannSolver

   abstract interface
      subroutine RiemannSolverFCN(QLeft, QRight, QbaseL, QbaseR, nHat, t1, t2, flux)
         use SMConstants
         use PhysicsStorage_CAA
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: QbaseL(1:NCONSB)
         real(kind=RP), intent(in)       :: QbaseR(1:NCONSB)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
      end subroutine RiemannSolverFCN
   end interface

   procedure(RiemannSolverFCN),      protected, pointer  :: RiemannSolver      => NULL()

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
         use RiemannSolvers_CAAKeywordsModule
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
!        ---------------------------------------------
!        Choose the Riemann solver (by default is LF)
!        ---------------------------------------------
         if (controlVariables % containsKey(RIEMANN_SOLVER_NAME_KEY)) then

            keyword = controlVariables % stringValueForKey(RIEMANN_SOLVER_NAME_KEY, KEYWORD_LENGTH)
            call toLower(keyword)

            select case (keyword)
            case (RIEMANN_ROE_NAME)
               RiemannSolver => RoeRiemannSolver
               whichRiemannSolver = RIEMANN_ROE

            case (RIEMANN_LXF_NAME)
               RiemannSolver => LxFRiemannSolver
               whichRiemannSolver = RIEMANN_LXF

            ! case (RIEMANN_RUSANOV_NAME)
            !   RiemannSolver => RusanovRiemannSolver
            !    whichRiemannSolver = RIEMANN_RUSANOV

            case (RIEMANN_CENTRAL_NAME)
               RiemannSolver => CentralRiemannSolver
               whichRiemannSolver = RIEMANN_CENTRAL

            case default
               print*, "Riemann Solver not recognized."
               errorMessage(STD_OUT)
               error stop

            end select

         else
!
!           Select LF by default
!           ---------------------
            RiemannSolver => LxFRiemannSolver
            whichRiemannSolver = RIEMANN_LXF

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
      end subroutine SetRiemannSolver

      subroutine DescribeRiemannSolver
!
!        -------
!        Modules
!        -------
         use RiemannSolvers_CAAKeywordsModule

         select case (whichRiemannSolver)
         case (RIEMANN_ROE)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Roe"

         case (RIEMANN_LXF)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Lax-Friedrichs"

         ! case (RIEMANN_RUSANOV)
            ! write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Rusanov"

         case (RIEMANN_CENTRAL)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Central"

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
      subroutine CentralRiemannSolver(QLeft, QRight, QbaseL, QbaseR, nHat, t1, t2, flux)
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: QbaseL(1:NCONSB)
         real(kind=RP), intent(in)       :: QbaseR(1:NCONSB)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: rhoL, uL, vL, wL, pL, a2L
         real(kind=RP) :: rhoR, uR, vR, wR, pR, a2R
         real(kind=RP) :: QLRot(NCONS), QRRot(NCONS) 
         real(kind=RP) :: QBaseLRot(NCONSB), QbaseRRot(NCONSB) 
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         rhoL = QLeft(1)
         rhoR = QRight(1)

         uL = QLeft(2) * nHat(1) + QLeft(3) * nHat(2) + QLeft(4) * nHat(3)
         vL = QLeft(2) * t1(1)   + QLeft(3) * t1(2)   + QLeft(4) * t1(3)
         wL = QLeft(2) * t2(1)   + QLeft(3) * t2(2)   + QLeft(4) * t2(3)

         uR = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)
         vR = QRight(2) * t1(1)   + QRight(3) * t1(2)   + QRight(4) * t1(3)
         wR = QRight(2) * t2(1)   + QRight(3) * t2(2)   + QRight(4) * t2(3)

         pL = QLeft(5)
         pR = QRight(5)
!
         QLRot = (/ rhoL, uL, vL, wL, pL /)
         QRRot = (/ rhoR, uR, vR, wR, pR /)
!
         rhoL = QbaseL(1)
         rhoR = QbaseR(1)

         uL = QbaseL(2) * nHat(1) + QbaseL(3) * nHat(2) + QbaseL(4) * nHat(3)
         vL = QbaseL(2) * t1(1)   + QbaseL(3) * t1(2)   + QbaseL(4) * t1(3)
         wL = QbaseL(2) * t2(1)   + QbaseL(3) * t2(2)   + QbaseL(4) * t2(3)

         uR = QbaseR(2) * nHat(1) + QbaseR(3) * nHat(2) + QbaseR(4) * nHat(3)
         vR = QbaseR(2) * t1(1)   + QbaseR(3) * t1(2)   + QbaseR(4) * t1(3)
         wR = QbaseR(2) * t2(1)   + QbaseR(3) * t2(2)   + QbaseR(4) * t2(3)

         pL = QbaseL(5)
         pR = QbaseR(5)

         a2L = QbaseL(6)
         a2R = QbaseR(6)
!
         QbaseLRot = (/ rhoL, uL, vL, wL, pL, a2L /)
         QbaseRRot = (/ rhoR, uR, vR, wR, pR, a2R /)
!        Perform the average using the averaging function
!        ------------------------------------------------
         call StandardAverage(QLRot, QRRot, QBaseLRot, QbaseRRot, flux)
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

      end subroutine CentralRiemannSolver
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LxFRiemannSolver( QLeft, QRight, QbaseL, QbaseR, nHat, t1, t2, flux )
         implicit none
!
!        ---------
!        Arguments
!        ---------
!
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: QbaseL(1:NCONSB)
         real(kind=RP), intent(in)       :: QbaseR(1:NCONSB)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local Variables
!        ---------------
!
         real(kind=RP) :: rhoL, uL, vL, wL, pL, aL
         real(kind=RP) :: rhoR, uR, vR, wR, pR, aR
         real(kind=RP) :: QLRot(NCONS), QRRot(NCONS) 
         real(kind=RP) :: QBaseLRot(NCONSB), QbaseRRot(NCONSB) 
         real(kind=RP)  :: lambda, stab(NCONS)
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         rhoL = QbaseL(1)
         rhoR = QbaseR(1)

         uL = QbaseL(2) * nHat(1) + QbaseL(3) * nHat(2) + QbaseL(4) * nHat(3)
         vL = QbaseL(2) * t1(1)   + QbaseL(3) * t1(2)   + QbaseL(4) * t1(3)
         wL = QbaseL(2) * t2(1)   + QbaseL(3) * t2(2)   + QbaseL(4) * t2(3)

         uR = QbaseR(2) * nHat(1) + QbaseR(3) * nHat(2) + QbaseR(4) * nHat(3)
         vR = QbaseR(2) * t1(1)   + QbaseR(3) * t1(2)   + QbaseR(4) * t1(3)
         wR = QbaseR(2) * t2(1)   + QbaseR(3) * t2(2)   + QbaseR(4) * t2(3)

         pL = QbaseL(5)
         pR = QbaseR(5)
!
         QbaseLRot = (/ rhoL, uL, vL, wL, pL, QBaseL(6) /)
         QbaseRRot = (/ rhoR, uR, vR, wR, pR, QBaseR(6) /)
!
         ! speed of sound of base flow
         aL = sqrt(QbaseL(6))
         aR = sqrt(QbaseR(6))
!
         rhoL = QLeft(1)
         rhoR = QRight(1)

         uL = QLeft(2) * nHat(1) + QLeft(3) * nHat(2) + QLeft(4) * nHat(3)
         vL = QLeft(2) * t1(1)   + QLeft(3) * t1(2)   + QLeft(4) * t1(3)
         wL = QLeft(2) * t2(1)   + QLeft(3) * t2(2)   + QLeft(4) * t2(3)

         uR = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)
         vR = QRight(2) * t1(1)   + QRight(3) * t1(2)   + QRight(4) * t1(3)
         wR = QRight(2) * t2(1)   + QRight(3) * t2(2)   + QRight(4) * t2(3)

         pL = QLeft(5)
         pR = QRight(5)
!
         QLRot = (/ rhoL, uL, vL, wL, pL /)
         QRRot = (/ rhoR, uR, vR, wR, pR /)
!
!        Eigenvalues: lambda = max(|uL|,|uR|) + max(aL,aR)
!        -----------
         ! lambda = max(abs(rhouL*invRhoL) + aL,abs(rhouR*invRhoR) + aR) ! this was the NS version
         lambda = max(abs(uL), abs(uR)) + max(aL,aR)
!
!        ****************
!        Compute the flux
!        ****************
!
!        Perform the average using the averaging function
!        ------------------------------------------------
         call StandardAverage(QLRot, QRRot, QBaseLRot, QbaseRRot, flux)
!
!        Compute the Lax-Friedrichs stabilization
!        ----------------------------------------
         stab = 0.5_RP * lambda * (QRRot - QLRot)
!
!        Compute the flux: apply the lambda stabilization here.
!        ----------------
         flux = flux - lambdaStab * stab
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

      END SUBROUTINE LxFRiemannSolver
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE RoeRiemannSolver( QLeft, QRight, QbaseL, QbaseR, nHat, t1, t2, flux )
!
!        **************************************************************
!           This Roe Riemann solver implementation does not uses the
!           averaging function, nor the lambda stabilization
!        **************************************************************
!
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: QbaseL(1:NCONSB)
         real(kind=RP), intent(in)       :: QbaseR(1:NCONSB)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
         ! REAL(KIND=RP) :: rho , rhou , rhov , rhow  , rhoe
         ! REAL(KIND=RP) :: rhon, rhoun, rhovn, rhown , rhoen
         ! REAL(KIND=RP) :: ul  , vl   , wl   , pleft , ql  , hl  , betal
         ! REAL(KIND=RP) :: ur  , vr   , wr   , pright, qr  , hr  , betar
         ! REAL(KIND=RP) :: rtd , utd  , vtd  , wtd   , htd , atd2, atd, qtd
         ! REAL(KIND=RP) :: dw1 , sp1  , sp1m , hd1m  , eta1, udw1, rql
         ! REAL(KIND=RP) :: dw4 , sp4  , sp4p , hd4   , eta4, udw4, rqr
         ! REAL(KIND=RP)                   :: ds = 1.0_RP

         !associate ( gamma => thermodynamics % gamma )

         !rho  = QLeft(1)
         !rhou = QLeft(2)
         !rhov = QLeft(3)
         !rhow = QLeft(4)
         !rhoe = QLeft(5)

         !rhon  = QRight(1)
         !rhoun = QRight(2)
         !rhovn = QRight(3)
         !rhown = QRight(4)
         !rhoen = QRight(5)

         !ul = rhou/rho
         !vl = rhov/rho
         !wl = rhow/rho
         !pleft = (gamma-1._RP)*(rhoe - 0.5_RP/rho*                        &
        !&                           (rhou**2 + rhov**2 + rhow**2 ))
!!
         !ur = rhoun/rhon
         !vr = rhovn/rhon
         !wr = rhown/rhon
         !pright = (gamma-1._RP)*(rhoen - 0.5_RP/rhon*                    &
        !&                           (rhoun**2 + rhovn**2+ rhown**2))
!!
         !ql = nHat(1)*ul + nHat(2)*vl + nHat(3)*wl
         !qr = nHat(1)*ur + nHat(2)*vr + nHat(3)*wr
         !hl = 0.5_RP*(ul*ul + vl*vl + wl*wl) +                               &
        !&                 gamma/(gamma-1._RP)*pleft/rho
         !hr = 0.5_RP*(ur*ur + vr*vr + wr*wr) +                               &
        !&                  gamma/(gamma-1._RP)*pright/rhon
!!
!!        ---------------------
!!        Square root averaging
!!        ---------------------
!!
         !rtd = sqrt(rho*rhon)
         !betal = rho/(rho + rtd)
         !betar = 1._RP - betal
         !utd = betal*ul + betar*ur
         !vtd = betal*vl + betar*vr
         !wtd = betal*wl + betar*wr
         !htd = betal*hl + betar*hr
         !atd2 = (gamma-1._RP)*(htd - 0.5_RP*(utd*utd + vtd*vtd + wtd*wtd))
         !atd = sqrt(atd2)
         !qtd = utd*nHat(1) + vtd*nHat(2)  + wtd*nHat(3)
!!
         !IF(qtd >= 0.0_RP)     THEN

         !   dw1 = 0.5_RP*((pright - pleft)/atd2 - (qr - ql)*rtd/atd)
         !   sp1 = qtd - atd
         !   sp1m = min(sp1,0.0_RP)
         !   hd1m = ((gamma+1._RP)/4._RP*atd/rtd)*dw1
         !   eta1 = max(-abs(sp1) - hd1m,0.0_RP)
         !   udw1 = dw1*(sp1m - 0.5_RP*eta1)
         !   rql = rho*ql
         !   flux(1) = ds*(rql + udw1)
         !   flux(2) = ds*(rql*ul + pleft*nHat(1) + udw1*(utd - atd*nHat(1)))
         !   flux(3) = ds*(rql*vl + pleft*nHat(2) + udw1*(vtd - atd*nHat(2)))
         !   flux(4) = ds*(rql*wl + pleft*nHat(3) + udw1*(wtd - atd*nHat(3)))
         !   flux(5) = ds*(rql*hl + udw1*(htd - qtd*atd))

         !ELSE

         !   dw4 = 0.5_RP*((pright - pleft)/atd2 + (qr - ql)*rtd/atd)
         !   sp4 = qtd + atd
         !   sp4p = max(sp4,0.0_RP)
         !   hd4 = ((gamma+1._RP)/4._RP*atd/rtd)*dw4
         !   eta4 = max(-abs(sp4) + hd4,0.0_RP)
         !   udw4 = dw4*(sp4p + 0.5_RP*eta4)
         !   rqr = rhon*qr
         !   flux(1) = ds*(rqr - udw4)
         !   flux(2) = ds*(rqr*ur + pright*nHat(1) - udw4*(utd + atd*nHat(1)))
         !   flux(3) = ds*(rqr*vr + pright*nHat(2) - udw4*(vtd + atd*nHat(2)))
         !   flux(4) = ds*(rqr*wr + pright*nHat(3) - udw4*(wtd + atd*nHat(3)))
         !   flux(5) = ds*(rqr*hr - udw4*(htd + qtd*atd))
         !ENDIF

         !end associate
print*, "Roe Riemann solver not implemented"
error stop

      END SUBROUTINE RoeRiemannSolver

!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      !SUBROUTINE RusanovRiemannSolver( QLeft, QRight, nHat, t1, t2, flux )

      !   IMPLICIT NONE
!!
!!        ---------
!!        Arguments
!!        ---------
!!
      !   real(kind=RP), intent(in)       :: QLeft(1:NCONS)
      !   real(kind=RP), intent(in)       :: QRight(1:NCONS)
      !   real(kind=RP), intent(in)       :: nHat(1:NDIM)
      !   real(kind=RP), intent(in)       :: t1(1:NDIM)
      !   real(kind=RP), intent(in)       :: t2(1:NDIM)
      !   real(kind=RP), intent(out)      :: flux(1:NCONS)
!!
!!        ---------------
!!        Local Variables
!!        ---------------
!!
!!
      !   REAL(KIND=RP) :: rho , rhou , rhov , rhow  , rhoe
      !   REAL(KIND=RP) :: rhon, rhoun, rhovn, rhown , rhoen
      !   REAL(KIND=RP) :: ul  , vl   , wl   , pleft , ql  , hl  , betal, al, al2
      !   REAL(KIND=RP) :: ur  , vr   , wr   , pright, qr  , hr  , betar, ar, ar2
      !   REAL(KIND=RP) :: rtd , utd  , vtd  , wtd   , htd , atd2, atd, qtd
      !   REAL(KIND=RP) :: dw1 , sp1  , sp1m , hd1m  , eta1, udw1, rql
      !   REAL(KIND=RP) :: dw4 , sp4  , sp4p , hd4   , eta4, udw4, rqr
      !   REAL(KIND=RP)                   :: ds = 1.0_RP

      !   REAL(KIND=RP) :: smax, smaxL, smaxR
      !   REAL(KIND=RP) :: Leigen(2), Reigen(2)

      !   associate ( gamma => thermodynamics % gamma )

      !   rho  = QLeft(1)
      !   rhou = QLeft(2)
      !   rhov = QLeft(3)
      !   rhow = QLeft(4)
      !   rhoe = QLeft(5)

      !   rhon  = QRight(1)
      !   rhoun = QRight(2)
      !   rhovn = QRight(3)
      !   rhown = QRight(4)
      !   rhoen = QRight(5)

      !   ul = rhou/rho
      !   vl = rhov/rho
      !   wl = rhow/rho
      !   pleft = (gamma-1._RP)*(rhoe - 0.5_RP/rho*                        &
      !  &                           (rhou**2 + rhov**2 + rhow**2 ))
!!
      !   ur = rhoun/rhon
      !   vr = rhovn/rhon
      !   wr = rhown/rhon
      !   pright = (gamma-1._RP)*(rhoen - 0.5_RP/rhon*                    &
      !  &                           (rhoun**2 + rhovn**2+ rhown**2))
!!
      !   ql = nHat(1)*ul + nHat(2)*vl + nHat(3)*wl
      !   qr = nHat(1)*ur + nHat(2)*vr + nHat(3)*wr
      !   hl = 0.5_RP*(ul*ul + vl*vl + wl*wl) +                               &
      !  &                 gamma/(gamma-1._RP)*pleft/rho
      !   hr = 0.5_RP*(ur*ur + vr*vr + wr*wr) +                               &
      !  &                  gamma/(gamma-1._RP)*pright/rhon
!!
!!        ---------------------
!!        Square root averaging
!!        ---------------------
!!
      !   rtd = sqrt(rho*rhon)
      !   betal = rho/(rho + rtd)
      !   betar = 1._RP - betal
      !   utd = betal*ul + betar*ur
      !   vtd = betal*vl + betar*vr
      !   wtd = betal*wl + betar*wr
      !   htd = betal*hl + betar*hr
      !   atd2 = (gamma-1._RP)*(htd - 0.5_RP*(utd*utd + vtd*vtd + wtd*wtd))
      !   atd = sqrt(atd2)
      !   qtd = utd*nHat(1) + vtd*nHat(2)  + wtd*nHat(3)
      !   !Rusanov
      !   ar2 = (gamma-1.d0)*(hr - 0.5d0*(ur*ur + vr*vr + wr*wr))
      !   al2 = (gamma-1.d0)*(hl - 0.5d0*(ul*ul + vl*vl + wl*wl))
      !   ar = sqrt(ar2)
      !   al = sqrt(al2)
!!
      !   rql = rho*ql
      !   rqr = rhon*qr
      !   flux(1) = ds*(rql + rqr)
      !   flux(2) = ds*(rql*ul + pleft*nHat(1) + rqr*ur + pright*nHat(1))
      !   flux(3) = ds*(rql*vl + pleft*nHat(2) + rqr*vr + pright*nHat(2))
      !   flux(4) = ds*(rql*wl + pleft*nHat(3) + rqr*wr + pright*nHat(3))
      !   flux(5) = ds*(rql*hl + rqr*hr)

      !   smax = MAX(ar+ABS(qr),al+ABS(ql))

      !   flux = (flux - ds*smax*(QRight-QLeft))/2.d0

      !   RETURN

      !   end associate

      !END SUBROUTINE RusanovRiemannSolver


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
!
!////////////////////////////////////////////////////////////////////////////////////////////
      subroutine StandardAverage(QLeft, QRight, QbaseL, QbaseR, flux)
!
!        *********************************************************************
!           Computes the standard average of the two states:
!              F* = {{F}} = 0.5 * (FL + FR)
!
!           State vectors are rotated.
!        *********************************************************************
!
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: QbaseL(1:NCONSB)
         real(kind=RP), intent(in)       :: QbaseR(1:NCONSB)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP), dimension(1:NCONS)     :: fL, fR

!
!        Compute the flux
!        ----------------
         call APEFlux1D(QRight, fR, Qbase=QbaseR)
         call APEFlux1D(QLeft,  fL, Qbase=QbaseL)

         flux = 0.5_RP*(fR+fL)

      end subroutine StandardAverage

end module RiemannSolvers_CAA
