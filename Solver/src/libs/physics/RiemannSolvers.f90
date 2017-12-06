!
!//////////////////////////////////////////////////////
!
!   @File:    RiemannSolvers.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Wed Dec  6 17:42:26 2017
!   @Last revision date: Wed Dec  6 17:59:35 2017
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 8cc17733385cd74924da75c3c1fd7cd89dd71fd4
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module RiemannSolvers
   use SMConstants
   use PhysicsStorage

   private 
   public RiemannSolver, SetRiemannSolver

   abstract interface
      subroutine RiemannSolverFCN(QLeft, QRight, nHat, flux)
         use SMConstants
         use PhysicsStorage
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
      end subroutine RiemannSolverFCN
   end interface

   procedure(RiemannSolverFCN), pointer   :: RiemannSolver => NULL()

   contains
      SUBROUTINE SetRiemannSolver(which)
         IMPLICIT NONE
         integer, intent(in) :: which

         select case ( which )
         case ( ROE )
            RiemannSolver => RoeSolver
         case (LXF)
            RiemannSolver => LxFSolver
         case (RUSANOV)
            RiemannSolver => RusanovSolver
         case (DUCROS)
            RiemannSolver => DucrosSolver
         case (MORINISHI)
            RiemannSolver => MorinishiSolver
         case (PIROZZOLI)
            RiemannSolver => PirozzoliSolver
         case (KENNEDYGRUBER)
            RiemannSolver => KennedyGruberSolver
         case default
            PRINT *, "Undefined choice of Riemann Solver. Abort"
            errorMessage(STD_OUT)
            STOP
         end select
      
      END SUBROUTINE SetRiemannSolver
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE RoeSolver( QLeft, QRight, nHat, flux )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
         REAL(KIND=RP) :: rho , rhou , rhov , rhow  , rhoe
         REAL(KIND=RP) :: rhon, rhoun, rhovn, rhown , rhoen
         REAL(KIND=RP) :: ul  , vl   , wl   , pleft , ql  , hl  , betal
         REAL(KIND=RP) :: ur  , vr   , wr   , pright, qr  , hr  , betar
         REAL(KIND=RP) :: rtd , utd  , vtd  , wtd   , htd , atd2, atd, qtd
         REAL(KIND=RP) :: dw1 , sp1  , sp1m , hd1m  , eta1, udw1, rql
         REAL(KIND=RP) :: dw4 , sp4  , sp4p , hd4   , eta4, udw4, rqr
         REAL(KIND=RP)                   :: ds = 1.0_RP
      
         associate ( gamma => thermodynamics % gamma )
            
         rho  = QLeft(1)
         rhou = QLeft(2)
         rhov = QLeft(3)
         rhow = QLeft(4)
         rhoe = QLeft(5)
   
         rhon  = QRight(1)
         rhoun = QRight(2)
         rhovn = QRight(3)
         rhown = QRight(4)
         rhoen = QRight(5)
         
         ul = rhou/rho 
         vl = rhov/rho 
         wl = rhow/rho 
         pleft = (gamma-1._RP)*(rhoe - 0.5_RP/rho*                        &
        &                           (rhou**2 + rhov**2 + rhow**2 )) 
!
         ur = rhoun/rhon 
         vr = rhovn/rhon 
         wr = rhown/rhon 
         pright = (gamma-1._RP)*(rhoen - 0.5_RP/rhon*                    &
        &                           (rhoun**2 + rhovn**2+ rhown**2)) 
!
         ql = nHat(1)*ul + nHat(2)*vl + nHat(3)*wl
         qr = nHat(1)*ur + nHat(2)*vr + nHat(3)*wr
         hl = 0.5_RP*(ul*ul + vl*vl + wl*wl) +                               &
        &                 gamma/(gamma-1._RP)*pleft/rho 
         hr = 0.5_RP*(ur*ur + vr*vr + wr*wr) +                               &
        &                  gamma/(gamma-1._RP)*pright/rhon 
!
!        ---------------------
!        Square root averaging  
!        ---------------------
!
         rtd = sqrt(rho*rhon) 
         betal = rho/(rho + rtd) 
         betar = 1._RP - betal 
         utd = betal*ul + betar*ur 
         vtd = betal*vl + betar*vr 
         wtd = betal*wl + betar*wr 
         htd = betal*hl + betar*hr 
         atd2 = (gamma-1._RP)*(htd - 0.5_RP*(utd*utd + vtd*vtd + wtd*wtd)) 
         atd = sqrt(atd2) 
         qtd = utd*nHat(1) + vtd*nHat(2)  + wtd*nHat(3)
!
         IF(qtd >= 0.0_RP)     THEN
   
            dw1 = 0.5_RP*((pright - pleft)/atd2 - (qr - ql)*rtd/atd) 
            sp1 = qtd - atd 
            sp1m = min(sp1,0.0_RP) 
            hd1m = ((gamma+1._RP)/4._RP*atd/rtd)*dw1 
            eta1 = max(-abs(sp1) - hd1m,0.0_RP) 
            udw1 = dw1*(sp1m - 0.5_RP*eta1) 
            rql = rho*ql 
            flux(1) = ds*(rql + udw1) 
            flux(2) = ds*(rql*ul + pleft*nHat(1) + udw1*(utd - atd*nHat(1))) 
            flux(3) = ds*(rql*vl + pleft*nHat(2) + udw1*(vtd - atd*nHat(2))) 
            flux(4) = ds*(rql*wl + pleft*nHat(3) + udw1*(wtd - atd*nHat(3))) 
            flux(5) = ds*(rql*hl + udw1*(htd - qtd*atd)) 
   
         ELSE 
   
            dw4 = 0.5_RP*((pright - pleft)/atd2 + (qr - ql)*rtd/atd) 
            sp4 = qtd + atd 
            sp4p = max(sp4,0.0_RP) 
            hd4 = ((gamma+1._RP)/4._RP*atd/rtd)*dw4 
            eta4 = max(-abs(sp4) + hd4,0.0_RP) 
            udw4 = dw4*(sp4p + 0.5_RP*eta4) 
            rqr = rhon*qr 
            flux(1) = ds*(rqr - udw4) 
            flux(2) = ds*(rqr*ur + pright*nHat(1) - udw4*(utd + atd*nHat(1))) 
            flux(3) = ds*(rqr*vr + pright*nHat(2) - udw4*(vtd + atd*nHat(2))) 
            flux(4) = ds*(rqr*wr + pright*nHat(3) - udw4*(wtd + atd*nHat(3))) 
            flux(5) = ds*(rqr*hr - udw4*(htd + qtd*atd)) 
         ENDIF

         end associate
         
      END SUBROUTINE RoeSolver
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LxFSolver( QLeft, QRight, nHat, flux ) 
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
         REAL(KIND=RP) :: rho , rhou , rhov  , rhow  , rhoe
         REAL(KIND=RP) :: rhon, rhoun, rhovn , rhown , rhoen
         REAL(KIND=RP) :: ul  , vl   , wl    , pleft , ql, cl
         REAL(KIND=RP) :: ur  , vr   , wr    , pright, qr, cr
         REAL(KIND=RP) :: sM
         REAL(KIND=RP), DIMENSION(N_EQN) :: FL, FR
         
         associate ( gamma => thermodynamics % gamma )
         
         rho  = QLeft(1)
         rhou = QLeft(2)
         rhov = QLeft(3)
         rhow = QLeft(4)
         rhoe = QLeft(5)
         
         rhon  = QRight(1)
         rhoun = QRight(2)
         rhovn = QRight(3)
         rhown = QRight(4)
         rhoen = QRight(5)
         
         ul = rhou/rho 
         vl = rhov/rho 
         wl = rhow/rho 
         pleft = (gamma-1.d0)*(rhoe - 0.5d0/rho*(rhou**2 + rhov**2 + rhow**2)) 
         
         ur = rhoun/rhon 
         vr = rhovn/rhon 
         wr = rhown/rhon 
         pright = (gamma-1.d0)*(rhoen - 0.5d0/rhon*(rhoun**2 + rhovn**2 + rhown**2)) 
         
         ql = nHat(1)*ul + nHat(2)*vl + nHat(3)*wl 
         qr = nHat(1)*ur + nHat(2)*vr + nHat(3)*wr 
         cl = sqrt( gamma*pleft/rho )
         cr = sqrt( gamma*pright/rhon )
         
         FL(1) = rho*ql
         FL(2) = rhou*ql + pleft*nHat(1)
         FL(3) = rhov*ql + pleft*nHat(2)
         FL(4) = rhow*ql + pleft*nHat(3)
         FL(5) = (rhoe + pleft)*ql
         
         FR(1) = rhon*qr
         FR(2) = rhoun*qr + pright*nHat(1)
         FR(3) = rhovn*qr + pright*nHat(2)
         FR(4) = rhown*qr + pright*nHat(3)
         FR(5) = (rhoen + pright)*qr
         
         sM = MAX( ABS(ql) + cl, ABS(qr) + cr )
         
         flux = 0.5_RP * ( FL + FR - (sM+lambdaStab)*(QRight - QLeft) )      
         
         end associate
      
      END SUBROUTINE LxFSolver
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE RusanovSolver( QLeft, QRight, nHat, flux )
      
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
         REAL(KIND=RP) :: rho , rhou , rhov , rhow  , rhoe
         REAL(KIND=RP) :: rhon, rhoun, rhovn, rhown , rhoen
         REAL(KIND=RP) :: ul  , vl   , wl   , pleft , ql  , hl  , betal, al, al2
         REAL(KIND=RP) :: ur  , vr   , wr   , pright, qr  , hr  , betar, ar, ar2
         REAL(KIND=RP) :: rtd , utd  , vtd  , wtd   , htd , atd2, atd, qtd
         REAL(KIND=RP) :: dw1 , sp1  , sp1m , hd1m  , eta1, udw1, rql
         REAL(KIND=RP) :: dw4 , sp4  , sp4p , hd4   , eta4, udw4, rqr
         REAL(KIND=RP)                   :: ds = 1.0_RP
         
         REAL(KIND=RP) :: smax, smaxL, smaxR
         REAL(KIND=RP) :: Leigen(2), Reigen(2)
      
         associate ( gamma => thermodynamics % gamma ) 

         rho  = QLeft(1)
         rhou = QLeft(2)
         rhov = QLeft(3)
         rhow = QLeft(4)
         rhoe = QLeft(5)
   
         rhon  = QRight(1)
         rhoun = QRight(2)
         rhovn = QRight(3)
         rhown = QRight(4)
         rhoen = QRight(5)
   
         ul = rhou/rho 
         vl = rhov/rho 
         wl = rhow/rho 
         pleft = (gamma-1._RP)*(rhoe - 0.5_RP/rho*                        &
        &                           (rhou**2 + rhov**2 + rhow**2 )) 
!
         ur = rhoun/rhon 
         vr = rhovn/rhon 
         wr = rhown/rhon 
         pright = (gamma-1._RP)*(rhoen - 0.5_RP/rhon*                    &
        &                           (rhoun**2 + rhovn**2+ rhown**2)) 
!
         ql = nHat(1)*ul + nHat(2)*vl + nHat(3)*wl
         qr = nHat(1)*ur + nHat(2)*vr + nHat(3)*wr
         hl = 0.5_RP*(ul*ul + vl*vl + wl*wl) +                               &
        &                 gamma/(gamma-1._RP)*pleft/rho 
         hr = 0.5_RP*(ur*ur + vr*vr + wr*wr) +                               &
        &                  gamma/(gamma-1._RP)*pright/rhon 
!
!        ---------------------
!        Square root averaging  
!        ---------------------
!
         rtd = sqrt(rho*rhon) 
         betal = rho/(rho + rtd) 
         betar = 1._RP - betal 
         utd = betal*ul + betar*ur 
         vtd = betal*vl + betar*vr 
         wtd = betal*wl + betar*wr 
         htd = betal*hl + betar*hr 
         atd2 = (gamma-1._RP)*(htd - 0.5_RP*(utd*utd + vtd*vtd + wtd*wtd)) 
         atd = sqrt(atd2) 
         qtd = utd*nHat(1) + vtd*nHat(2)  + wtd*nHat(3)
         !Rusanov
         ar2 = (gamma-1.d0)*(hr - 0.5d0*(ur*ur + vr*vr + wr*wr)) 
         al2 = (gamma-1.d0)*(hl - 0.5d0*(ul*ul + vl*vl + wl*wl)) 
         ar = sqrt(ar2)
         al = sqrt(al2)
!           
         rql = rho*ql 
         rqr = rhon*qr             
         flux(1) = ds*(rql + rqr) 
         flux(2) = ds*(rql*ul + pleft*nHat(1) + rqr*ur + pright*nHat(1)) 
         flux(3) = ds*(rql*vl + pleft*nHat(2) + rqr*vr + pright*nHat(2))
         flux(4) = ds*(rql*wl + pleft*nHat(3) + rqr*wr + pright*nHat(3)) 
         flux(5) = ds*(rql*hl + rqr*hr) 

         smax = MAX(ar+ABS(qr),al+ABS(ql))

         flux = (flux - ds*smax*(QRight-QLeft))/2.d0

         RETURN 

         end associate
         
      END SUBROUTINE RusanovSolver           

      subroutine DucrosSolver(QLeft, QRight, nHat, flux) 
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, uL, vL, wL, pL, aL, unL
         real(kind=RP)     :: invRhoR, uR, vR, wR, pR, aR, unR
         real(kind=RP)     :: f(NCONS), g(NCONS), h(NCONS)

         invRhoL = 1.0_RP / QLeft(IRHO)   ; invRhoR = 1.0_RP / QRight(IRHO)
         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)
   
         pL = thermodynamics % GammaMinus1 * ( QLeft(IRHOE) - 0.5_RP * (   QLeft(IRHOU) * uL &
                                                                      + QLeft(IRHOV) * vL &
                                                                      + QLeft(IRHOW) * wL ))

         pR = thermodynamics % GammaMinus1 * ( QRight(IRHOE) - 0.5_RP * (   QRight(IRHOU) * uR &
                                                                      + QRight(IRHOV) * vR &
                                                                      + QRight(IRHOW) * wR ))
         aL = sqrt(thermodynamics % gamma * pL * invRhoL)
         aR = sqrt(thermodynamics % gamma * pR * invRhoR)
         unL = sum( (/uL,vL,wL/) * nHat ) ; unR = sum( (/uR,vR,wR/) * nHat )
!
!        Compute the flux
!        ----------------
         f(IRHO)  = 0.25_RP * ( QLeft(IRHO) + QRight(IRHO) ) * (uL + uR)
         f(IRHOU) = 0.25_RP * ( QLeft(IRHOU) + QRight(IRHOU) ) * (uL + uR) + 0.5_RP * (pL + pR)
         f(IRHOV) = 0.25_RP * ( QLeft(IRHOV) + QRight(IRHOV) ) * (uL + uR)
         f(IRHOW) = 0.25_RP * ( QLeft(IRHOW) + QRight(IRHOW) ) * (uL + uR)
         f(IRHOE) = 0.25_RP * ( QLeft(IRHOE) + pL + QRight(IRHOE) + pR ) * (uL + uR)

         g(IRHO)  = 0.25_RP * ( QLeft(IRHO) + QRight(IRHO) ) * (vL + vR)
         g(IRHOU) = 0.25_RP * ( QLeft(IRHOU) + QRight(IRHOU) ) * (vL + vR)
         g(IRHOV) = 0.25_RP * ( QLeft(IRHOV) + QRight(IRHOV) ) * (vL + vR) + 0.5_RP * (pL + pR)
         g(IRHOW) = 0.25_RP * ( QLeft(IRHOW) + QRight(IRHOW) ) * (vL + vR)
         g(IRHOE) = 0.25_RP * ( QLeft(IRHOE) + pL + QRight(IRHOE) + pR ) * (vL + vR)

         h(IRHO)  = 0.25_RP * ( QLeft(IRHO) + QRight(IRHO) ) * (wL + wR)
         h(IRHOU) = 0.25_RP * ( QLeft(IRHOU) + QRight(IRHOU) ) * (wL + wR)
         h(IRHOV) = 0.25_RP * ( QLeft(IRHOV) + QRight(IRHOV) ) * (wL + wR)
         h(IRHOW) = 0.25_RP * ( QLeft(IRHOW) + QRight(IRHOW) ) * (wL + wR) + 0.5_RP * (pL + pR)
         h(IRHOE) = 0.25_RP * ( QLeft(IRHOE) + pL + QRight(IRHOE) + pR ) * (wL + wR)
!
!        Compute the sharp flux
!        ----------------------         
         flux = f*nHat(IX) + g*nHat(IY) + h*nHat(IZ) - 0.5_RP * lambdaStab * max(abs(unL)+aL,abs(unR)+aR) * (QRight-QLeft)

      end subroutine DucrosSolver

      subroutine MorinishiSolver(QLeft,QRight,nHat,flux) 
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(NDIM)
         real(kind=RP), intent(out)      :: flux(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, uL, vL, wL, pL, hL, aL, unL
         real(kind=RP)     :: invRhoR, uR, vR, wR, pR, hR, aR, unR
         real(kind=RP)     :: f(NCONS), g(NCONS), h(NCONS)

         invRhoL = 1.0_RP / QLeft(IRHO)   ; invRhoR = 1.0_RP / QRight(IRHO)
         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)
   
         pL = thermodynamics % GammaMinus1 * ( QLeft(IRHOE) - 0.5_RP * (   QLeft(IRHOU) * uL &
                                                                      + QLeft(IRHOV) * vL &
                                                                      + QLeft(IRHOW) * wL ))

         pR = thermodynamics % GammaMinus1 * ( QRight(IRHOE) - 0.5_RP * (   QRight(IRHOU) * uR &
                                                                      + QRight(IRHOV) * vR &
                                                                      + QRight(IRHOW) * wR ))
         aL = sqrt(thermodynamics % gamma * pL * invRhoL)
         aR = sqrt(thermodynamics % gamma * pR * invRhoR)
         unL = sum( (/uL,vL,wL/) * nHat ) ; unR = sum( (/uR,vR,wR/) * nHat )
!
!        Here the enthalpy does not contain the kinetic energy
!        -----------------------------------------------------
         hL = dimensionless % cp * pL  ; hR = dimensionless % cp * pR 
!
!        Compute the flux
!        ----------------
         f(IRHO)  = 0.5_RP * ( QLeft(IRHOU) + QRight(IRHOU) )
         f(IRHOU) = 0.25_RP * ( QLeft(IRHOU) + QRight(IRHOU) ) * ( uL + uR ) + 0.5_RP * ( pL + pR )
         f(IRHOV) = 0.25_RP * ( QLeft(IRHOU) + QRight(IRHOU) ) * ( vL + vR )
         f(IRHOW) = 0.25_RP * ( QLeft(IRHOU) + QRight(IRHOU) ) * ( wL + wR )
         f(IRHOE) = 0.5_RP * ( uL*hL + uR*hR) + 0.25_RP * ( QLeft(IRHOU)*uL + QRight(IRHOU)*uR ) * ( uL + uR ) &
                                              + 0.25_RP * ( QLeft(IRHOU)*vL + QRight(IRHOU)*vR ) * ( vL + vR ) &
                                              + 0.25_RP * ( QLeft(IRHOU)*wL + QRight(IRHOU)*wR ) * ( wL + wR ) &
                                              - 0.25_RP * ( QLeft(IRHOU)*POW2(uL) + QRight(IRHOU)*POW2(uR)   ) &
                                              - 0.25_RP * ( QLeft(IRHOU)*POW2(vL) + QRight(IRHOU)*POW2(vR)   ) &
                                              - 0.25_RP * ( QLeft(IRHOU)*POW2(wL) + QRight(IRHOU)*POW2(wR)   ) 

         g(IRHO)  = 0.5_RP * ( QLeft(IRHOV) + QRight(IRHOV) )
         g(IRHOU) = 0.25_RP * ( QLeft(IRHOV) + QRight(IRHOV) ) * ( uL + uR )
         g(IRHOV) = 0.25_RP * ( QLeft(IRHOV) + QRight(IRHOV) ) * ( vL + vR ) + 0.5_RP * ( pL + pR )
         g(IRHOW) = 0.25_RP * ( QLeft(IRHOV) + QRight(IRHOV) ) * ( wL + wR )
         g(IRHOE) = 0.5_RP * ( vL*hL + vR*hR) + 0.25_RP * ( QLeft(IRHOV)*uL + QRight(IRHOV)*uR ) * ( uL + uR ) &
                                              + 0.25_RP * ( QLeft(IRHOV)*vL + QRight(IRHOV)*vR ) * ( vL + vR ) &
                                              + 0.25_RP * ( QLeft(IRHOV)*wL + QRight(IRHOV)*wR ) * ( wL + wR ) &
                                              - 0.25_RP * ( QLeft(IRHOV)*POW2(uL) + QRight(IRHOV)*POW2(uR)   ) &
                                              - 0.25_RP * ( QLeft(IRHOV)*POW2(vL) + QRight(IRHOV)*POW2(vR)   ) &
                                              - 0.25_RP * ( QLeft(IRHOV)*POW2(wL) + QRight(IRHOV)*POW2(wR)   ) 

         h(IRHO)  = 0.5_RP * ( QLeft(IRHOW) + QRight(IRHOW) )
         h(IRHOU) = 0.25_RP * ( QLeft(IRHOW) + QRight(IRHOW) ) * ( uL + uR )
         h(IRHOV) = 0.25_RP * ( QLeft(IRHOW) + QRight(IRHOW) ) * ( vL + vR )
         h(IRHOW) = 0.25_RP * ( QLeft(IRHOW) + QRight(IRHOW) ) * ( wL + wR ) + 0.5_RP * ( pL + pR )
         h(IRHOE) = 0.5_RP * ( wL*hL + wR*hR) + 0.25_RP * ( QLeft(IRHOW)*uL + QRight(IRHOW)*uR ) * ( uL + uR ) &
                                              + 0.25_RP * ( QLeft(IRHOW)*vL + QRight(IRHOW)*vR ) * ( vL + vR ) &
                                              + 0.25_RP * ( QLeft(IRHOW)*wL + QRight(IRHOW)*wR ) * ( wL + wR ) &
                                              - 0.25_RP * ( QLeft(IRHOW)*POW2(uL) + QRight(IRHOW)*POW2(uR)   ) &
                                              - 0.25_RP * ( QLeft(IRHOW)*POW2(vL) + QRight(IRHOW)*POW2(vR)   ) &
                                              - 0.25_RP * ( QLeft(IRHOW)*POW2(wL) + QRight(IRHOW)*POW2(wR)   ) 

!
!        Compute the sharp flux
!        ----------------------         
         flux = f*nHat(IX) + g*nHat(IY) + h*nHat(IZ) - 0.5_RP * lambdaStab * max(abs(unL)+aL,abs(unR)+aR) * (QRight-QLeft)

      end subroutine MorinishiSolver

      subroutine KennedyGruberSolver(QLeft,QRight,nHat,flux) 
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(NDIM)
         real(kind=RP), intent(out)      :: flux(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, uL, vL, wL, pL, aL, unL
         real(kind=RP)     :: invRhoR, uR, vR, wR, pR, aR, unR
         real(kind=RP)     :: rho, u, v, w, e, p
         real(kind=RP)     :: f(NCONS), g(NCONS), h(NCONS)

         invRhoL = 1.0_RP / QLeft(IRHO)   ; invRhoR = 1.0_RP / QRight(IRHO)
         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)
   
         pL = thermodynamics % GammaMinus1 * ( QLeft(IRHOE) - 0.5_RP * (   QLeft(IRHOU) * uL &
                                                                      + QLeft(IRHOV) * vL &
                                                                      + QLeft(IRHOW) * wL ))

         pR = thermodynamics % GammaMinus1 * ( QRight(IRHOE) - 0.5_RP * (   QRight(IRHOU) * uR &
                                                                      + QRight(IRHOV) * vR &
                                                                      + QRight(IRHOW) * wR ))

         aL = sqrt(thermodynamics % gamma * pL * invRhoL)
         aR = sqrt(thermodynamics % gamma * pR * invRhoR)
         unL = sum( (/uL,vL,wL/) * nHat ) ; unR = sum( (/uR,vR,wR/) * nHat )

         rho = 0.5_RP * (QLeft(IRHO) + QRight(IRHO))
         u   = 0.5_RP * (uL + uR)
         v   = 0.5_RP * (vL + vR)
         w   = 0.5_RP * (wL + wR)
         p   = 0.5_RP * (pL + pR)
         e   = 0.5_RP * (QLeft(IRHOE)*invRhoL + QRight(IRHOE)*invRhoR)
!
!        Compute the flux
!        ----------------
         f(IRHO)  = rho * u
         f(IRHOU) = rho * u * u + p
         f(IRHOV) = rho * u * v
         f(IRHOW) = rho * u * w
         f(IRHOE) = rho * u * e + p * u
         
         g(IRHO)  = rho * v
         g(IRHOU) = rho * v * u
         g(IRHOV) = rho * v * v + p
         g(IRHOW) = rho * v * w
         g(IRHOE) = rho * v * e + p * v

         h(IRHO)  = rho * w
         h(IRHOU) = rho * w * u
         h(IRHOV) = rho * w * v
         h(IRHOW) = rho * w * w + p
         h(IRHOE) = rho * w * e + p * w
!
!        Compute the sharp flux
!        ----------------------         
         flux = f*nHat(IX) + g*nHat(IY) + h*nHat(IZ) - 0.5_RP * lambdaStab * max(abs(unL)+aL,abs(unR)+aR) * (QRight-QLeft)

      end subroutine KennedyGruberSolver

      subroutine PirozzoliSolver(QLeft,QRight,nHat,flux)
         use SMConstants
         use PhysicsStorage
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(NDIM)
         real(kind=RP), intent(out)      :: flux(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: invRhoL, uL, vL, wL, pL, aL, unL
         real(kind=RP)     :: invRhoR, uR, vR, wR, pR, aR, unR
         real(kind=RP)     :: rho, u, v, w, h, p
         real(kind=RP)     :: ff(NCONS), gg(NCONS), hh(NCONS)

         invRhoL = 1.0_RP / QLeft(IRHO)   ; invRhoR = 1.0_RP / QRight(IRHO)
         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)
   
         pL = thermodynamics % GammaMinus1 * ( QLeft(IRHOE) - 0.5_RP * (   QLeft(IRHOU) * uL &
                                                                      + QLeft(IRHOV) * vL &
                                                                      + QLeft(IRHOW) * wL ))

         pR = thermodynamics % GammaMinus1 * ( QRight(IRHOE) - 0.5_RP * (   QRight(IRHOU) * uR &
                                                                      + QRight(IRHOV) * vR &
                                                                      + QRight(IRHOW) * wR ))
         aL = sqrt(thermodynamics % gamma * pL * invRhoL)
         aR = sqrt(thermodynamics % gamma * pR * invRhoR)
         unL = sum( (/uL,vL,wL/) * nHat ) ; unR = sum( (/uR,vR,wR/) * nHat )

         rho = 0.5_RP * (QLeft(IRHO) + QRight(IRHO))
         u   = 0.5_RP * (uL + uR)
         v   = 0.5_RP * (vL + vR)
         w   = 0.5_RP * (wL + wR)
         p   = 0.5_RP * (pL + pR)
         h   = 0.5_RP * ((QLeft(IRHOE)+pL)*invRhoL + (QRight(IRHOE)+pR)*invRhoR)
!
!        Compute the flux
!        ----------------
         ff(IRHO)  = rho * u
         ff(IRHOU) = rho * u * u + p
         ff(IRHOV) = rho * u * v
         ff(IRHOW) = rho * u * w
         ff(IRHOE) = rho * u * h
         
         gg(IRHO)  = rho * v
         gg(IRHOU) = rho * v * u
         gg(IRHOV) = rho * v * v + p
         gg(IRHOW) = rho * v * w
         gg(IRHOE) = rho * v * h

         hh(IRHO)  = rho * w
         hh(IRHOU) = rho * w * u
         hh(IRHOV) = rho * w * v
         hh(IRHOW) = rho * w * w + p
         hh(IRHOE) = rho * w * h 
!
!        Compute the sharp flux
!        ----------------------         
         flux = ff*nHat(IX) + gg*nHat(IY) + hh*nHat(IZ) - 0.5_RP * lambdaStab * max(abs(unL)+aL,abs(unR)+aR) * (QRight-QLeft)

      end subroutine PirozzoliSolver
end module RiemannSolvers
