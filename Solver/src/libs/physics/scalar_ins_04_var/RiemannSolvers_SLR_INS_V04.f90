#include "Includes.h"

! module RiemannSolvers_SLR_INS_V04KeywordsModule


! end module RiemannSolvers_SLR_INS_V04KeywordsModule

module RiemannSolvers_SLR_INS_V04


    use SMConstants
    use Physics_SLR_INS_V04
    use PhysicsStorage_SLR_INS_V04
    use FluidData_SLR_INS_V04

    implicit none

    private

    ! public whichAverage, whichRiemannSolver
    ! public SetRiemannSolver, DescribeRiemannSolver
    ! public TwoPointFlux, RiemannSolver_dFdQ
    public LxFRiemannSolver, AveragedStates
    public LxFRiemannSolverTest, LxFRiemannSolver_BC_Test
    public LxFRiemannSolver_burges
 
    real(RP)           :: lambdaStab = 1.0_RP
 !
 !  ========
    contains
 !  ========
 !
    !
!////////////////////////////////////////////////////////////////////////
!
    SUBROUTINE LxFRiemannSolver( QLeft, QRight, nHat, t1, t2, flux )
        implicit none
!
!        ---------
!        Arguments
!        ---------
!
        real(kind=RP), intent(in)       :: QLeft(1:N_INS)
        real(kind=RP), intent(in)       :: QRight(1:N_INS)
        real(kind=RP), intent(in)       :: nHat(1:NDIM)
        real(kind=RP), intent(in)       :: t1(1:NDIM)
        real(kind=RP), intent(in)       :: t2(1:NDIM)
        real(kind=RP), intent(out)      :: flux(1:N_INS)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
        real(kind=RP)  :: uxL, uyL, uzL
        real(kind=RP)  :: uxR, uyR, uzR
        real(kind=RP)  :: QLRot(N_INS), QRRot(N_INS)
        real(kind=RP)  :: invRhoL, invRhoR
        real(kind=RP)  :: lambda, stab(N_INS)
        real(kind=RP)  :: uLMax, uRMax

! !        Rotate the variables to the face local frame using normal and tangent vectors
! !        -----------------------------------------------------------------------------

        uxL = QLeft (1) * nHat(1) + QLeft (2) * nHat(2) + QLeft (3) * nHat(3)
        uxR = QRight(1) * nHat(1) + QRight(2) * nHat(2) + QRight(3) * nHat(3)

        uyL = QLeft (1) * t1(1) + QLeft (2) * t1(2) + QLeft (3) * t1(3)
        uyR = QRight(1) * t1(1) + QRight(2) * t1(2) + QRight(3) * t1(3)

        uzL = QLeft (1) * t2(1) + QLeft (2) * t2(2) + QLeft (3) * t2(3)
        uzR = QRight(1) * t2(1) + QRight(2) * t2(2) + QRight(3) * t2(3)

! !     Eigenvalues: lambda = max(|uL|,|uR|) + max(aL,aR)
        ! uLMax = abs(uxL)
        ! uRMax = abs(uxR)
        uLMax = abs(uxL) + abs(uyL) + abs(uzL)
        uRMax = abs(uxR) + abs(uyR) + abs(uzR)
! !        -----------
        lambda = max(uLMax, uRMax)
! !
!        ****************
!        Compute the flux
!        ****************
!
!        Perform the average using the averaging function
!        ------------------------------------------------
        QLRot = (/ uxL, uyL, uzL /)
        QRRot = (/ uxR, uyR, uzR /)
        call AveragedStates(QLeft=QLRot, QRight=QRRot, flux = flux)
! !
! !        Compute the Lax-Friedrichs stabilization
! !        ----------------------------------------
        stab = 0.5_RP * lambda * (QRRot - QLRot)
! !
! !        Compute the flux: apply the lambda stabilization here.
! !        ----------------
        ! flux = flux 
        flux = flux - lambdaStab * stab
! !
! !        ************************************************
! !        Return momentum equations to the cartesian frame
! !        ************************************************
! !
        flux(1:N_INS) = nHat*flux(1) + t1*flux(2) + t2*flux(3)

        ! end associate

     END SUBROUTINE LxFRiemannSolver
!////////////////////////////////////////////////////////////////////////
!
    SUBROUTINE LxFRiemannSolverTest_BK( QLeft, QRight, nHat, t1, t2, flux )
        implicit none
!
!        ---------
!        Arguments
!        ---------
!
        real(kind=RP), intent(in)       :: QLeft(1:N_INS)
        real(kind=RP), intent(in)       :: QRight(1:N_INS)
        real(kind=RP), intent(in)       :: nHat(1:NDIM)
        real(kind=RP), intent(in)       :: t1(1:NDIM)
        real(kind=RP), intent(in)       :: t2(1:NDIM)
        real(kind=RP), intent(out)      :: flux(1:N_INS)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
        real(kind=RP)  :: uxL, uyL, uzL
        real(kind=RP)  :: uxR, uyR, uzR

        real(kind=RP)  :: fxL(1:N_INS), fyL(1:N_INS), fzL(1:N_INS)
        real(kind=RP)  :: fxR(1:N_INS), fyR(1:N_INS), fzR(1:N_INS)
        real(kind=RP)  :: fxL_NR, fxR_NR
        real(kind=RP)  :: fyL_NR, fyR_NR
        real(kind=RP)  :: fzL_NR, fzR_NR

        real(kind=RP)  :: QLRot(N_INS), QRRot(N_INS)
        real(kind=RP)  :: invRhoL, invRhoR
        real(kind=RP)  :: lambda, stab(N_INS)
        real(kind=RP)  :: uLMax, uRMax

        real(kind=RP)  :: lambdaX, lambdaY, lambdaZ

        real(kind=RP) :: f(1:N_INS), g(1:N_INS), h(1:N_INS)

        
!         associate(gamma => thermodynamics % gamma, gm1 => thermodynamics % gammaMinus1)
! !
! !        Rotate the variables to the face local frame using normal and tangent vectors
! !        -----------------------------------------------------------------------------
!         rhoL = QLeft(1)                  ; rhoR = QRight(1)
!         invRhoL = 1.0_RP/ rhoL           ; invRhoR = 1.0_RP / rhoR


        lambdaX = 0.4 * nHat(1) + 0.0 * nHat(2) + 0.8 * nHat(3) 

        f(1) = lambdaX * ( QLeft(1) + QRight(1) ) * 0.5_RP
        f(2) = lambdaX * ( QLeft(2) + QRight(2) ) * 0.5_RP
        f(3) = lambdaX * ( QLeft(3) + QRight(3) ) * 0.5_RP
        g(1) = lambdaX * ( QLeft(1) + QRight(1) ) * 0.5_RP
        g(2) = lambdaX * ( QLeft(2) + QRight(2) ) * 0.5_RP
        g(3) = lambdaX * ( QLeft(3) + QRight(3) ) * 0.5_RP
        h(1) = lambdaX * ( QLeft(1) + QRight(1) ) * 0.5_RP
        h(2) = lambdaX * ( QLeft(2) + QRight(2) ) * 0.5_RP
        h(3) = lambdaX * ( QLeft(3) + QRight(3) ) * 0.5_RP

        lambdaX = 0.4 * nHat(1) + 0.0 * nHat(2) + 0.8 * nHat(3) 
        ! lambdaY = 0.4 * nHat(1) + 0.0 * nHat(2) + 0.8 * nHat(3) 
        ! lambdaZ = 0.4 * nHat(1) + 0.0 * nHat(2) + 0.8 * nHat(3) 

        ! write (*,*) "lambda X Y Z = ", lambdaX,lambdaY, lambdaZ



        ! f(1) = f(1) - 0.5_RP * lambdaX * ( QRight(1) - QLeft(1)  )
        ! f(2) = f(2) - 0.5_RP * lambdaY * ( QRight(2) - QLeft(2)  )
        ! f(3) = f(3) - 0.5_RP * lambdaZ * ( QRight(3) - QLeft(3)  )
        ! g(1) = g(1) - 0.5_RP * lambdaX * ( QRight(1) - QLeft(1)  )
        ! g(2) = g(2) - 0.5_RP * lambdaY * ( QRight(2) - QLeft(2)  )
        ! g(3) = g(3) - 0.5_RP * lambdaZ * ( QRight(3) - QLeft(3)  )
        ! h(1) = h(1) - 0.5_RP * lambdaX * ( QRight(1) - QLeft(1)  )
        ! h(2) = h(2) - 0.5_RP * lambdaY * ( QRight(2) - QLeft(2)  )
        ! h(3) = h(3) - 0.5_RP * lambdaZ * ( QRight(3) - QLeft(3)  )



        flux = f * nHat(1) + g*nHat(2) + h*nHat(3)

     END SUBROUTINE LxFRiemannSolverTest_BK


    SUBROUTINE LxFRiemannSolverTest( QLeft, QRight, nHat, t1, t2, flux )
        implicit none
!
!        ---------
!        Arguments
!        ---------
!
        real(kind=RP), intent(in)       :: QLeft(1:N_INS)
        real(kind=RP), intent(in)       :: QRight(1:N_INS)
        real(kind=RP), intent(in)       :: nHat(1:NDIM)
        real(kind=RP), intent(in)       :: t1(1:NDIM)
        real(kind=RP), intent(in)       :: t2(1:NDIM)
        real(kind=RP), intent(out)      :: flux(1:N_INS)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
        real(kind=RP)  :: uxL, uyL, uzL
        real(kind=RP)  :: uxR, uyR, uzR

        real(kind=RP)  :: fxL(1:N_INS), fyL(1:N_INS), fzL(1:N_INS)
        real(kind=RP)  :: fxR(1:N_INS), fyR(1:N_INS), fzR(1:N_INS)
        real(kind=RP)  :: fxL_NR, fxR_NR
        real(kind=RP)  :: fyL_NR, fyR_NR
        real(kind=RP)  :: fzL_NR, fzR_NR

        real(kind=RP)  :: QLRot(N_INS), QRRot(N_INS)
        real(kind=RP)  :: invRhoL, invRhoR
        real(kind=RP)  :: lambda, stab(N_INS)
        real(kind=RP)  :: uLMax, uRMax

        real(kind=RP)  :: lambdaX, lambdaY, lambdaZ

        real(kind=RP) :: f(1:N_INS), g(1:N_INS), h(1:N_INS)


        real(kind=RP) :: cx, cy, cz, cc


        cx = 0.5_RP
        cy = 0.0_RP
        cz = 0.5_RP


        lambdaX = cx * nHat(1) + cy * nHat(2) + cz * nHat(3) 
        lambdaY = cx * nHat(1) + cy * nHat(2) + cz * nHat(3) 
        lambdaZ = cx * nHat(1) + cy * nHat(2) + cz * nHat(3) 
        ! lambdaX = cx

        f(1) = lambdaX * ( QLeft(1) + QRight(1) ) * 0.5_RP
        ! f(2) = lambdaX * ( QLeft(1) + QRight(1) ) * 0.5_RP
        ! f(3) = lambdaX * ( QLeft(1) + QRight(1) ) * 0.5_RP
        f(2) = lambdaX * ( QLeft(2) + QRight(2) ) * 0.5_RP
        f(3) = lambdaX * ( QLeft(3) + QRight(3) ) * 0.5_RP
        g(1) = lambdaY * ( QLeft(1) + QRight(1) ) * 0.5_RP
        g(2) = lambdaY * ( QLeft(2) + QRight(2) ) * 0.5_RP
        g(3) = lambdaY * ( QLeft(3) + QRight(3) ) * 0.5_RP
        h(1) = lambdaZ * ( QLeft(1) + QRight(1) ) * 0.5_RP
        h(2) = lambdaZ * ( QLeft(2) + QRight(2) ) * 0.5_RP
        h(3) = lambdaZ * ( QLeft(3) + QRight(3) ) * 0.5_RP

        ! lambdaX = cx * nHat(1) + cy * nHat(2) + cz * nHat(3) 
        ! lambdaY = cx * nHat(1) + cy * nHat(2) + cz * nHat(3) 
        ! lambdaZ = cx * nHat(1) + cy * nHat(2) + cz * nHat(3) 



        f(1) = f(1) - 0.5_RP * abs(lambdaX) * ( QRight(1) - QLeft(1)  )
        f(2) = f(2) - 0.5_RP * abs(lambdaX) * ( QRight(2) - QLeft(2)  )
        f(3) = f(3) - 0.5_RP * abs(lambdaX) * ( QRight(3) - QLeft(3)  )
        g(1) = g(1) - 0.5_RP * abs(lambdaY) * ( QRight(1) - QLeft(1)  )
        g(2) = g(2) - 0.5_RP * abs(lambdaY) * ( QRight(2) - QLeft(2)  )
        g(3) = g(3) - 0.5_RP * abs(lambdaY) * ( QRight(3) - QLeft(3)  )
        h(1) = h(1) - 0.5_RP * abs(lambdaZ) * ( QRight(1) - QLeft(1)  )
        h(2) = h(2) - 0.5_RP * abs(lambdaZ) * ( QRight(2) - QLeft(2)  )
        h(3) = h(3) - 0.5_RP * abs(lambdaZ) * ( QRight(3) - QLeft(3)  )

        ! write (*,*) "f================",f
        ! write (*,*) "g================",g
        ! write (*,*) "h================",h
        ! write (*,*) "nHat=============",nHat
        ! write (*,*) "QLeft =============",QLeft
        ! write (*,*) "QRight=============",QRight
        ! write (*,*) "Qleft - QRight=============",QLeft-QRight



        flux = f * nHat(1) + g*nHat(2) + h*nHat(3)
        ! flux = f
        ! write (*,*) "f * nHat(1) + g*nHat(2) + h*nHat(3)",f * nHat(1) + g*nHat(2) + h*nHat(3)
        ! write (*,*) "flux =============",flux

     END SUBROUTINE LxFRiemannSolverTest

    SUBROUTINE LxFRiemannSolver_burges( QLeft, QRight, nHat, t1, t2, flux )
        implicit none
!
!        ---------
!        Arguments
!        ---------
!
        real(kind=RP), intent(in)       :: QLeft(1:N_INS)
        real(kind=RP), intent(in)       :: QRight(1:N_INS)
        real(kind=RP), intent(in)       :: nHat(1:NDIM)
        real(kind=RP), intent(in)       :: t1(1:NDIM)
        real(kind=RP), intent(in)       :: t2(1:NDIM)
        real(kind=RP), intent(out)      :: flux(1:N_INS)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
        real(kind=RP)  :: un_L, ut1L, ut2L
        real(kind=RP)  :: un_R, ut1R, ut2R

        real(kind=RP)  :: fxL(1:N_INS), fyL(1:N_INS), fzL(1:N_INS)
        real(kind=RP)  :: fxR(1:N_INS), fyR(1:N_INS), fzR(1:N_INS)
        real(kind=RP)  :: fxL_NR, fxR_NR
        real(kind=RP)  :: fyL_NR, fyR_NR
        real(kind=RP)  :: fzL_NR, fzR_NR

        real(kind=RP)  :: QLRot(N_INS), QRRot(N_INS)
        real(kind=RP)  :: invRhoL, invRhoR
        real(kind=RP)  :: lambda, stab(N_INS)
        real(kind=RP)  :: uLMax, uRMax

        real(kind=RP)  :: lambdaX, lambdaY, lambdaZ

        real(kind=RP) :: f(1:N_INS), g(1:N_INS), h(1:N_INS)


        real(kind=RP) :: cx, cy, cz, cc

        un_L = QLeft(1)  * nHat(1) + QLeft(2)  * nHat(2) + QLeft(3)  * nHat(3)
        ut1L = QLeft(1)  * t1  (1) + QLeft(2)  * t1  (2) + QLeft(3)  * t1  (3)
        ut2L = QLeft(1)  * t2  (1) + QLeft(2)  * t2  (2) + QLeft(3)  * t2  (3)


        un_R = QRight(1) * nHat(1) + QRight(2) * nHat(2) + QRight(3) * nHat(3)
        ut1R = QRight(1) * t1  (1) + QRight(2) * t1  (2) + QRight(3) * t1  (3)
        ut2R = QRight(1) * t2  (1) + QRight(2) * t2  (2) + QRight(3) * t2  (3)


        f(1) = 0.5_RP * ( un_L * un_L + un_R * un_R )
        f(2) = 0.5_RP * ( un_L * ut1L + un_R * ut1R )
        f(3) = 0.5_RP * ( un_L * ut2L + un_R * ut2R )

        lambda = max(abs(un_L), abs(un_R))

        ! 0.5_RP * lambda * (QRRot - QLRot)

        f(1) = f(1) 
        f(2) = f(2) 
        f(3) = f(3) 
        ! f(1) = f(1) - 0.5_RP * lambda * (un_R - un_L)
        ! f(2) = f(2) - 0.5_RP * lambda * (ut1R - ut1L)
        ! f(3) = f(3) - 0.5_RP * lambda * (ut2R - ut2L)
        ! write (*,*) "lambda       ******* =============",lambda
        ! write (*,*) "QLeft    ************ =============",QLeft
        ! write (*,*) "QRight   ************ =============",QRight
        ! write (*,*) "f    ************ =============",f
        ! write (*,*) "nhat ************ =============",nHat

        ! flux = f * nHat + 0.5_RP * (f(2) * t1 + f(3) * t2) - 0.5_RP * lambda * (QRight - QLeft)

        flux = f(1) * nHat + f(2) * t1 + f(3) * t2
        ! flux = f
        ! write (*,*) "f * nHat(1) + g*nHat(2) + h*nHat(3)",f * nHat(1) + g*nHat(2) + h*nHat(3)
        ! write (*,*) "flux =============",flux

     END SUBROUTINE LxFRiemannSolver_burges



!
     SUBROUTINE LxFRiemannSolver_BC_Test( QLeft, QRight, nHat, t1__, t2__, flux )
        implicit none
!
!        ---------
!        Arguments
!        ---------
!
        real(kind=RP), intent(in)       :: QLeft(1:N_INS)
        real(kind=RP), intent(in)       :: QRight(1:N_INS)
        real(kind=RP), intent(in)       :: nHat(1:NDIM)
        real(kind=RP), intent(in)       :: t1__(1:NDIM)
        real(kind=RP), intent(in)       :: t2__(1:NDIM)
        real(kind=RP), intent(out)      :: flux(1:N_INS)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
        real(kind=RP)  :: uxL, uyL, uzL
        real(kind=RP)  :: uxR, uyR, uzR

        real(kind=RP)  :: fxL(1:N_INS), fyL(1:N_INS), fzL(1:N_INS)
        real(kind=RP)  :: fxR(1:N_INS), fyR(1:N_INS), fzR(1:N_INS)
        real(kind=RP)  :: fxL_NR, fxR_NR
        real(kind=RP)  :: fyL_NR, fyR_NR
        real(kind=RP)  :: fzL_NR, fzR_NR

        real(kind=RP)  :: QLRot(N_INS), QRRot(N_INS)
        real(kind=RP)  :: invRhoL, invRhoR
        real(kind=RP)  :: lambda, stab(N_INS)
        real(kind=RP)  :: uLMax, uRMax

        real(kind=RP)  :: lambdaX, lambdaY, lambdaZ

        real(kind=RP) :: f(1:N_INS), g(1:N_INS), h(1:N_INS)

        real(kind=RP) :: cx, cy, cz, cc


        cx = 0.4_RP
        cy = 0.0_RP
        cz = 0.0_RP

        cc =  cx * nHat(1) + cy * nHat(2) + cz * nHat(3)
        uxL =  QLeft(1) * nHat(1) + QLeft(2) * nHat(2) + QLeft(3) * nHat(3)
        uyL =  QLeft(1) * t1__(1) + QLeft(2) * t1__(2) + QLeft(3) * t1__(3)
        uzL =  QLeft(1) * t2__(1) + QLeft(2) * t2__(2) + QLeft(3) * t2__(3)
        
        uxR =  QRight(1) * nHat(1) + QRight(2) * nHat(2) + QRight(3) * nHat(3)
        uyR =  QRight(1) * t1__(1) + QRight(2) * t1__(2) + QRight(3) * t1__(3)
        uzR =  QRight(1) * t2__(1) + QRight(2) * t2__(2) + QRight(3) * t2__(3)

        write (*,*) " cc,cx,cy,cz, =========== ",cc, cx,cy,cz
        write (*,*) " QLeft   ,  =========== ", QLeft
        write (*,*) " QRight  ,  =========== ", QRight
        write (*,*) " uxL,uyL,uzL , =========== ",uxL,uyL,uzL 
        write (*,*) " uxR,uyR,uzR , =========== ",uxR,uyR,uzR 
        write (*,*) " nHat , =========== ",nHat
        write (*,*) " t1__ , =========== ",t1__
        write (*,*) " t2__ , =========== ",t2__


        if ( cc .ge. 0.0_RP ) then
                f(1) = cx * QLeft(1)
                f(2) = cx * QLeft(2)
                f(3) = cx * QLeft(3)
                g(1) = cy * QLeft(1)
                g(2) = cy * QLeft(2)
                g(3) = cy * QLeft(3)
                h(1) = cz * QLeft(1)
                h(2) = cz * QLeft(2)
                h(3) = cz * QLeft(3)
        else 
                f(1) = cx * QRight(1)
                f(2) = cx * QRight(2)
                f(3) = cx * QRight(3)
                g(1) = cy * QRight(1)
                g(2) = cy * QRight(2)
                g(3) = cy * QRight(3)
                h(1) = cz * QRight(1)
                h(2) = cz * QRight(2)
                h(3) = cz * QRight(3)
        endif
  
        ! flux = f * nHat(1) + g*nHat(2) + h*nHat(3)
        flux = f * nHat(1) + g*nHat(2) + h*nHat(3)
        ! flux = f
        ! flux = 0

     END SUBROUTINE LxFRiemannSolver_BC_Test



     !
!      subroutine AveragedStates(QLeft, QRight, f, g, h)
subroutine AveragedStates(QLeft, QRight, flux)
        !
        !        *********************************************************************
        !           Computes the standard average of the two states:
        !              F* = {{F}} = 0.5 * (FL + FR)
        !
        !           State vectors are rotated.
        !        *********************************************************************
        !
                 implicit none
                 real(kind=RP), intent(in)       :: QLeft(1:N_INS)
                 real(kind=RP), intent(in)       :: QRight(1:N_INS)
                 real(kind=RP), intent(out)      :: flux(1:N_INS)
                !  real(kind=RP), intent(out)      :: f(1:N_INS), g(1:N_INS), h(1:N_INS)
        !
        !
        !        ---------------
        !        Local variables
        !        ---------------
        !
                 real(kind=RP)     :: uL, vL, wL
                 real(kind=RP)     :: uR, vR, wR
        
                 uL = QLeft(1)      ; uR = QRight(1)
                 vL = QLeft(2)      ; vR = QRight(2)
                 wL = QLeft(3)      ; wR = QRight(3)
        ! ! !
        ! !        Compute the flux
        ! ! !        ----------------
        ! !          flux(IRHO)  = 0.5_RP * ( QLeft(IRHOU) + QRight(IRHOU) )
                 flux(1) = 0.5_RP * ( QLeft(1) * uL + QRight(1) * uR )
                 flux(2) = 0.5_RP * ( QLeft(1) * vL + QRight(1) * vR )
                 flux(3) = 0.5_RP * ( QLeft(1) * wL + QRight(1) * wR )
        !          flux(IRHOE) = 0.5_RP * ( uL*(QLeft(IRHOE) + pL) + uR*(QRight(IRHOE) + pR) )
        
    end subroutine AveragedStates


    subroutine SkewSymmetric2Average(QLeft, QRight, f, g, h)
        !
        !        *********************************************************************
        !        *********************************************************************
        !
                 use Physics_iNS, only: iEulerXFlux
                 implicit none
                 real(kind=RP), intent(in)       :: QLeft(1:N_INS)
                 real(kind=RP), intent(in)       :: QRight(1:N_INS)
                 real(kind=RP), intent(out)      :: f(1:N_INS), g(1:N_INS), h(1:N_INS)
        !
        !        ---------------
        !        Local variables
        !        ---------------
        !
                 real(kind=RP)  :: invRhoL, invRhoR
                 real(kind=RP)  :: rho, u, v, w, p
        
        
                 u   = 0.5_RP * (QLeft(1) + QRight(1))
                 v   = 0.5_RP * (QLeft(2) + QRight(2))
                 w   = 0.5_RP * (QLeft(3) + QRight(3))
        
                 f(1) = f(1)*u
                 f(2) = f(1)*v
                 f(3) = f(1)*w
        
                 g(1) = g(1)*u
                 g(2) = g(1)*v
                 g(3) = g(1)*w
        
                 h(1) = h(1)*u
                 h(2) = h(1)*v
                 h(3) = h(1)*w
        
        end subroutine SkewSymmetric2Average


end module RiemannSolvers_SLR_INS_V04