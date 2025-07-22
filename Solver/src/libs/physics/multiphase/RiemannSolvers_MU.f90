#include "Includes.h"
module RiemannSolvers_MUKeywordsModule

     integer,                       parameter :: KEYWORD_LENGTH          = 132
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_SOLVER_NAME_KEY = "riemann solver"
!
!    --------------------------
!    Riemann solver definitions
!    --------------------------
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_CENTRAL_NAME = "central"
     character(len=KEYWORD_LENGTH), parameter :: RIEMANN_EXACT_NAME   = "exact"

     enum, bind(C)
        enumerator :: RIEMANN_CENTRAL = 1, RIEMANN_EXACT
     end enum
     !$acc declare copyin(RIEMANN_CENTRAL,RIEMANN_EXACT)

end module RiemannSolvers_MUKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
!
module RiemannSolvers_MU
   use SMConstants
   use Physics_MU
   use PhysicsStorage_MU
   use VariableConversion_MU
   use FluidData_MU

   implicit none

   private
   public whichRiemannSolver
   public SetRiemannSolver, DescribeRiemannSolver
   public ExactRiemannSolver
   public RiemannSolver_Selector_MU

   abstract interface
   !    subroutine RiemannSolverFCN(QLeft, QRight, rhoL, rhoR, muL, muR, nHat, t1, t2, fL,fR)
   !       use SMConstants
   !       use PhysicsStorage_MU
   !       real(kind=RP), intent(in)       :: QLeft(1:NCONS)
   !       real(kind=RP), intent(in)       :: QRight(1:NCONS)
   !       real(kind=RP), intent(in)       :: rhoL
   !       real(kind=RP), intent(in)       :: rhoR
   !       real(kind=RP), intent(in)       :: muL
   !       real(kind=RP), intent(in)       :: muR
   !       real(kind=RP), intent(in)       :: nHat(1:NDIM)
   !       real(kind=RP), intent(in)       :: t1(1:NDIM)
   !       real(kind=RP), intent(in)       :: t2(1:NDIM)
   !       real(kind=RP), intent(out)      :: fL(1:NCONS)
   !       real(kind=RP), intent(out)      :: fR(1:NCONS)
   !    end subroutine RiemannSolverFCN
   end interface

   integer           :: whichRiemannSolver = -1
   !$acc declare copyin(whichRiemannSolver)
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
         use RiemannSolvers_MUKeywordsModule
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
               !RiemannSolver => CentralRiemannSolver
               whichRiemannSolver = RIEMANN_CENTRAL

            case(RIEMANN_EXACT_NAME)
               !RiemannSolver => ExactRiemannSolver
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
            !RiemannSolver => CentralRiemannSolver
            whichRiemannSolver = RIEMANN_EXACT

         end if

         !$acc update device(whichRiemannSolver)

      END SUBROUTINE SetRiemannSolver

      subroutine DescribeRiemannSolver
!
!        -------
!        Modules
!        -------
         use RiemannSolvers_MUKeywordsModule


         select case (whichRiemannSolver)
         case (RIEMANN_CENTRAL)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Central"

         case (RIEMANN_EXACT)
            write(STD_OUT,'(30X,A,A30,A)') "->","Riemann solver: ","Exact"

         end select

      end subroutine DescribeRiemannSolver

      subroutine RiemannSolver_Selector_MU(Nx, Ny, QLeft, QRight, rhoL, rhoR, muL, muR, nHat, t1, t2, fL, fR, invMa2L, invMa2R)
         !$acc routine vector
         use RiemannSolvers_MUKeywordsModule
         implicit none 
         integer , intent(in)            :: Nx, Ny
         real(kind=RP), intent(in)       :: QLeft(1:NCONS, 0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: QRight(1:NCONS, 0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: rhoL(0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: rhoR(0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: muL(0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: muR(0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: nHat(1:NDIM, 0:Nx, 0:Ny), t1(NDIM, 0:Nx, 0:Ny), t2(NDIM, 0:Nx, 0:Ny)
         real(kind=RP), intent(inout)    :: fL(1:NCONS, 0:Nx, 0:Ny)
         real(kind=RP), intent(inout)    :: fR(1:NCONS, 0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: invMa2L(0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: invMa2R(0:Nx, 0:Ny)


         select case (whichRiemannSolver)
         case (RIEMANN_CENTRAL)
            call CentralRiemannSolver(Nx, Ny, QLeft, QRight, rhoL, rhoR, muL, muR, nHat, t1, t2, fL, fR, invMa2L, invMa2R)
         case (RIEMANN_EXACT)
            call ExactRiemannSolver(Nx, Ny, QLeft, QRight, rhoL, rhoR, muL, muR, nHat, t1, t2, fL, fR, invMa2L, invMa2R)
         end select
         
      end subroutine RiemannSolver_Selector_MU
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
!        Riemann solvers
!        ---------------
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine CentralRiemannSolver(Nx, Ny, QLeft, QRight, rhoL, rhoR, muL, muR, nHat, t1, t2, fL, fR, invMa2L, invMa2R)
         !$acc routine vector
         implicit none
         integer , intent(in)            :: Nx, Ny
         real(kind=RP), intent(in)       :: QLeft(1:NCONS, 0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: QRight(1:NCONS, 0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: rhoL(0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: rhoR(0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: muL(0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: muR(0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: nHat(1:NDIM, 0:Nx, 0:Ny), t1(NDIM, 0:Nx, 0:Ny), t2(NDIM, 0:Nx, 0:Ny)
         real(kind=RP), intent(out)      :: fL(1:NCONS, 0:Nx, 0:Ny)
         real(kind=RP), intent(out)      :: fR(1:NCONS, 0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: invMa2L(0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: invMa2R(0:Nx, 0:Ny)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: cL, uL, vL, wL, pL, invSqrtRhoL
         real(kind=RP) :: cR, uR, vR, wR, pR, invSqrtRhoR
         integer :: i,j
         real(kind=RP) :: flux_rot_R(5), flux_rot_L(5)


         !$acc loop vector collapse(2) private(flux_rot_R, flux_rot_L)
         do j = 0, Ny 
         do i = 0, Nx 
!
!        Left state variables and fluxes
!        -------------------------------
            invSqrtRhoL = 1.0_RP / sqrt(rhoL(i,j))
            cL = QLeft(IMC,i,j)
            uL = invSqrtRhoL * (QLeft(IMSQRHOU,i,j) * nHat(1,i,j) + QLeft(IMSQRHOV,i,j) * nHat(2,i,j) + QLeft(IMSQRHOW,i,j) * nHat(3,i,j))
            vL = invSqrtRhoL * (QLeft(IMSQRHOU,i,j) * t1(1,i,j)   + QLeft(IMSQRHOV,i,j) * t1(2,i,j)   + QLeft(IMSQRHOW,i,j) * t1(3,i,j))
            wL = invSqrtRhoL * (QLeft(IMSQRHOU,i,j) * t2(1,i,j)   + QLeft(IMSQRHOV,i,j) * t2(2,i,j)   + QLeft(IMSQRHOW,i,j) * t2(3,i,j))
            pL = QLeft(IMP,i,j)

            flux_rot_L(IMC)      = uL*cL
            flux_rot_L(IMSQRHOU) = 0.5_RP*rhoL(i,j)*uL*uL + pL
            flux_rot_L(IMSQRHOV) = 0.5_RP*rhoL(i,j)*uL*vL
            flux_rot_L(IMSQRHOW) = 0.5_RP*rhoL(i,j)*uL*wL
            flux_rot_L(IMP)      = 0.0_RP

!
!        Right state variables and fluxes
!        --------------------------------
            invSqrtRhoR = 1.0_RP / sqrt(rhoR(i,j))
            cR = QRight(IMC,i,j)
            uR = invSqrtRhoR * (QRight(IMSQRHOU,i,j) * nHat(1,i,j) + QRight(IMSQRHOV,i,j) * nHat(2,i,j) + QRight(IMSQRHOW,i,j) * nHat(3,i,j))
            vR = invSqrtRhoR * (QRight(IMSQRHOU,i,j) * t1(1,i,j)   + QRight(IMSQRHOV,i,j) * t1(2,i,j)   + QRight(IMSQRHOW,i,j) * t1(3,i,j))
            wR = invSqrtRhoR * (QRight(IMSQRHOU,i,j) * t2(1,i,j)   + QRight(IMSQRHOV,i,j) * t2(2,i,j)   + QRight(IMSQRHOW,i,j) * t2(3,i,j))
            pR = QRight(IMP,i,j)

            flux_rot_R(IMC)      = uR*cR
            flux_rot_R(IMSQRHOU) = 0.5_RP*rhoR(i,j)*uR*uR + pR
            flux_rot_R(IMSQRHOV) = 0.5_RP*rhoR(i,j)*uR*vR
            flux_rot_R(IMSQRHOW) = 0.5_RP*rhoR(i,j)*uR*wR
            flux_rot_R(IMP)      = 0.0_RP
!
!        Perform the average and rotation
!        --------------------------------
            flux_rot_L(:) = 0.5_RP*(flux_rot_L(:) + flux_rot_R(:))
            flux_rot_R(:) = flux_rot_L(:)
!
!        Add the non-conservative term
!        -----------------------------
            flux_rot_L(IMSQRHOU) = flux_rot_L(IMSQRHOU) + 0.5_RP*cL*(muR(i,j)-muL(i,j)) + 0.25_RP*rhoL(i,j)*uL*(uR-uL)
            flux_rot_L(IMSQRHOV) = flux_rot_L(IMSQRHOV) + 0.25_RP*rhoL(i,j)*uL*(vR-vL)
            flux_rot_L(IMSQRHOW) = flux_rot_L(IMSQRHOW) + 0.25_RP*rhoL(i,j)*uL*(wR-wL)
            flux_rot_L(IMP)      = flux_rot_L(IMP)      + 0.5_RP*invMa2L(i,j)*(uR-uL)

            flux_rot_R(IMSQRHOU) = flux_rot_R(IMSQRHOU) + 0.5_RP*cR*(muL(i,j)-muR(i,j)) + 0.25_RP*rhoR(i,j)*uR*(uL-uR)
            flux_rot_R(IMSQRHOV) = flux_rot_R(IMSQRHOV) + 0.25_RP*rhoR(i,j)*uR*(vL-vR)
            flux_rot_R(IMSQRHOW) = flux_rot_R(IMSQRHOW) + 0.25_RP*rhoR(i,j)*uR*(wL-wR)
            flux_rot_R(IMP)      = flux_rot_R(IMP)      + 0.5_RP*invMa2R(i,j)*(uL-uR)

            
            fL(1,i,j) = flux_rot_L(1)
            fL(2,i,j) = nHat(1,i,j)*flux_rot_L(2) + t1(1,i,j)*flux_rot_L(3) + t2(1,i,j)*flux_rot_L(4)
            fL(3,i,j) = nHat(2,i,j)*flux_rot_L(2) + t1(2,i,j)*flux_rot_L(3) + t2(2,i,j)*flux_rot_L(4)
            fL(4,i,j) = nHat(3,i,j)*flux_rot_L(2) + t1(3,i,j)*flux_rot_L(3) + t2(3,i,j)*flux_rot_L(4)  
            fL(5,i,j) = flux_rot_L(5)

            fR(1,i,j) = flux_rot_R(1)
            fR(2,i,j) = nHat(1,i,j)*flux_rot_R(2) + t1(1,i,j)*flux_rot_R(3) + t2(1,i,j)*flux_rot_R(4)
            fR(3,i,j) = nHat(2,i,j)*flux_rot_R(2) + t1(2,i,j)*flux_rot_R(3) + t2(2,i,j)*flux_rot_R(4)
            fR(4,i,j) = nHat(3,i,j)*flux_rot_R(2) + t1(3,i,j)*flux_rot_R(3) + t2(3,i,j)*flux_rot_R(4)  
            fR(5,i,j) = flux_rot_R(5)
            
         enddo 
         enddo
         
      end subroutine CentralRiemannSolver

      subroutine ExactRiemannSolver(Nx, Ny, QLeft, QRight, rhoL, rhoR, muL, muR, nHat, t1, t2, fL, fR, invMa2L, invMa2R)
         !$acc routine vector
         implicit none
         integer , intent(in)            :: Nx, Ny
         real(kind=RP), intent(in)       :: QLeft(1:NCONS, 0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: QRight(1:NCONS, 0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: rhoL(0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: rhoR(0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: muL(0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: muR(0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: nHat(1:NDIM, 0:Nx, 0:Ny), t1(NDIM, 0:Nx, 0:Ny), t2(NDIM, 0:Nx, 0:Ny)
         real(kind=RP), intent(inout)    :: fL(1:NCONS, 0:Nx, 0:Ny)
         real(kind=RP), intent(inout)    :: fR(1:NCONS, 0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: invMa2L(0:Nx, 0:Ny)
         real(kind=RP), intent(in)       :: invMa2R(0:Nx, 0:Ny)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: cL,uL, vL, wL, pL, invRhoL, invSqrtRhoL, lambdaMinusL, lambdaPlusL
         real(kind=RP)  :: cR,uR, vR, wR, pR, invRhoR, invSqrtRhoR, lambdaMinusR, lambdaPlusR
         real(kind=RP)  :: rhoStarL, rhoStarR, uStar, pStar, rhoStar, vStar, wStar, cuStar, halfRhouStar
         real(kind=RP)  :: lambda_mu = 0.0_RP
         real(kind=RP)  :: flux_rot_R(5), flux_rot_L(5)

         integer :: i,j

         !$acc loop vector collapse(2) private(flux_rot_R, flux_rot_L)
         do j = 0,Ny
            do i = 0,Nx 
!  
!           Rotate the variables to the face local frame using normal and tangent vectors
!           -----------------------------------------------------------------------------
            invRhoL     = 1.0_RP / rhoL(i,j)
            invSqrtRhoL = sqrt(invRhoL)
            cL = QLeft(IMC,i,j)
            uL = invSqrtRhoL * (QLeft(IMSQRHOU,i,j) * nHat(1,i,j) + QLeft(IMSQRHOV,i,j) * nHat(2,i,j) + QLeft(IMSQRHOW,i,j) * nHat(3,i,j))
            vL = invSqrtRhoL * (QLeft(IMSQRHOU,i,j) * t1(1,i,j)   + QLeft(IMSQRHOV,i,j) * t1(2,i,j)   + QLeft(IMSQRHOW,i,j) * t1(3,i,j))
            wL = invSqrtRhoL * (QLeft(IMSQRHOU,i,j) * t2(1,i,j)   + QLeft(IMSQRHOV,i,j) * t2(2,i,j)   + QLeft(IMSQRHOW,i,j) * t2(3,i,j))
            pL = QLeft(IMP,i,j)

            invRhoR     = 1.0_RP / rhoR(i,j)
            invSqrtRhoR = sqrt(invRhoR)
            cR = QRight(IMC,i,j)
            uR = invSqrtRhoR * (QRight(IMSQRHOU,i,j) * nHat(1,i,j) + QRight(IMSQRHOV,i,j) * nHat(2,i,j) + QRight(IMSQRHOW,i,j) * nHat(3,i,j))
            vR = invSqrtRhoR * (QRight(IMSQRHOU,i,j) * t1(1,i,j)   + QRight(IMSQRHOV,i,j) * t1(2,i,j)   + QRight(IMSQRHOW,i,j) * t1(3,i,j))
            wR = invSqrtRhoR * (QRight(IMSQRHOU,i,j) * t2(1,i,j)   + QRight(IMSQRHOV,i,j) * t2(2,i,j)   + QRight(IMSQRHOW,i,j) * t2(3,i,j))
            pR = QRight(IMP,i,j)
!  
!           Compute the Star Region
!           -----------------------
            lambdaMinusR = 0.5_RP * (uR - sqrt(uR*uR + 4.0_RP*invMa2R(i,j)/rhoR(i,j)))
            lambdaPlusR  = 0.5_RP * (uR + sqrt(uR*uR + 4.0_RP*invMa2R(i,j)/rhoR(i,j)))

            lambdaMinusL = 0.5_RP * (uL - sqrt(uL*uL + 4.0_RP*invMa2L(i,j)/rhoL(i,j)))
            lambdaPlusL  = 0.5_RP * (uL + sqrt(uL*uL + 4.0_RP*invMa2L(i,j)/rhoL(i,j)))

            uStar = (pR-pL+rhoR(i,j)*uR*lambdaMinusR-rhoL(i,j)*uL*lambdaPlusL)/(rhoR(i,j)*lambdaMinusR - rhoL(i,j)*lambdaPlusL)
            pStar = pR + rhoR(i,j)*lambdaMinusR*(uR-uStar)
            rhoStarL = (rhoL(i,j)*lambdaPlusL)/(uStar-lambdaMinusL)
            rhoStarR = (rhoR(i,j)*lambdaMinusR)/(uStar - lambdaPlusR)

            if ( uStar .ge. 0.0_RP ) then
               rhoStar = rhoStarL
               vStar   = vL
               wStar   = wL

            else
               rhoStar = rhoStarR
               vStar   = vR
               wStar   = wR

            end if

            cuStar = 0.5_RP*(cL*uL + cR*uR)
            halfRhouStar = 0.5_RP*rhoStar*uStar
!
            ! Add first the common (conservative) part
            !fL = [cuStar+lambda_mu*(muL-muR), rhoStar*uStar*uStar + pStar, rhoStar*uStar*vStar, rhoStar*uStar*wStar, dimensionless % invMa2 * uStar]
            flux_rot_L = [cuStar+lambda_mu*(muL(i,j)-muR(i,j)), halfRhouStar*uStar + pStar, halfRhouStar*vStar, halfRhouStar*wStar, 0.0_RP] ! 0.5*(invMa2L+invMa2R) * uStar
            flux_rot_R = flux_rot_L
!  
!      -    Add the non--conservative part
            flux_rot_L = flux_rot_L + [0.0_RP, cL*0.5_RP*(muR(i,j)-muL(i,j))+ 0.5_RP*halfRhouStar*(uStar-uL),0.5_RP*halfRhouStar*(vStar-vL), 0.5_RP*halfRhouStar*(wStar-wL), (invMa2L(i,j))*(uStar-uL)]
            flux_rot_R = flux_rot_R + [0.0_RP, cR*0.5_RP*(muL(i,j)-muR(i,j))+ 0.5_RP*halfRhouStar*(uStar-uR),0.5_RP*halfRhouStar*(vStar-vR), 0.5_RP*halfRhouStar*(wStar-wR), (invMa2R(i,j))*(uStar-uR)]

!            ************************************************
!            Return momentum equations to the cartesian frame
!            ************************************************
!

            fL(1,i,j) = flux_rot_L(1)
            fL(2,i,j) = nHat(1,i,j)*flux_rot_L(2) + t1(1,i,j)*flux_rot_L(3) + t2(1,i,j)*flux_rot_L(4)
            fL(3,i,j) = nHat(2,i,j)*flux_rot_L(2) + t1(2,i,j)*flux_rot_L(3) + t2(2,i,j)*flux_rot_L(4)
            fL(4,i,j) = nHat(3,i,j)*flux_rot_L(2) + t1(3,i,j)*flux_rot_L(3) + t2(3,i,j)*flux_rot_L(4)  
            fL(5,i,j) = flux_rot_L(5)

            fR(1,i,j) = flux_rot_R(1)
            fR(2,i,j) = nHat(1,i,j)*flux_rot_R(2) + t1(1,i,j)*flux_rot_R(3) + t2(1,i,j)*flux_rot_R(4)
            fR(3,i,j) = nHat(2,i,j)*flux_rot_R(2) + t1(2,i,j)*flux_rot_R(3) + t2(2,i,j)*flux_rot_R(4)
            fR(4,i,j) = nHat(3,i,j)*flux_rot_R(2) + t1(3,i,j)*flux_rot_R(3) + t2(3,i,j)*flux_rot_R(4)  
            fR(5,i,j) = flux_rot_R(5)
            
            enddo 
         enddo

      end subroutine ExactRiemannSolver

end module RiemannSolvers_MU