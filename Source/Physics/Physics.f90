!      Physics.f90
!      Created: 2011-07-20 09:17:26 -0400 
!      By: David Kopriva
!      From DSEM Code
!
!!     The variable mappings for the Navier-Stokes Equations are
!!
!!              Q(1) = rho
!!              Q(2) = rhou
!!              Q(3) = rhov
!!              Q(4) = rhow
!!              Q(5) = rhoe
!!     Whereas the gradients are:
!!              grad(1) = grad(u)
!!              grad(2) = grad(v)
!!              grad(3) = grad(w)
!!              grad(4) = grad(T)
!
!////////////////////////////////////////////////////////////////////////
!    
!    ******
     MODULE PhysicsStorage
!    ******
!
     USE SMConstants
     
     IMPLICIT NONE
     SAVE
!
!    ----------------------------
!    Either NavierStokes or Euler
!    ----------------------------
!
     LOGICAL :: flowIsNavierStokes = .true.
!
!    --------------------------
!!   The sizes of the NS system
!    --------------------------
!
     INTEGER :: N_EQN = 5, N_GRAD_EQN = 4
!
!    ----------------------------------------
!!   The free-stream or reference mach number
!    ----------------------------------------
!
     REAL( KIND=RP ) :: mach 
!
!    ----------------------------------------
!!   The Reynolds number
!    ----------------------------------------
!
     REAL( KIND=RP ) :: RE 
!
!    ----------------------------------------
!!   The free-stream Angle of Attack
!    ----------------------------------------
!
     REAL( KIND=RP ) :: AOATheta, AOAPhi
!
!    ----------------------------------------
!!   The Prandtl number
!    ----------------------------------------
!
     REAL( KIND=RP ) :: PR 
!
!    ----------------------------------------
!!   The free-stream or reference temperature
!!   with default in R.
!    ----------------------------------------
!
     REAL( KIND=RP ) :: TRef 
!
!    --------------------------------------------
!!   The temperature scale in the Sutherland law:
!!   198.6 for temperatures in R, 110.3 for
!!   temperatures in K.
!    --------------------------------------------
!
     REAL( KIND=RP ) :: TScale
!
!    ------------------------------------------------
!!   The ratio of the scale and reference tempartures
!    ------------------------------------------------
!
     REAL( KIND=RP ) :: TRatio 
!
!    -------------
!!   The gas gamma
!    -------------
!
     REAL( KIND=RP ) :: gamma
!
!    ----------------------------------
!!   Other constants derived from gamma
!    ----------------------------------
!
     REAL( KIND=RP ) :: sqrtGamma          , gammaMinus1      , gammaMinus1Div2
     REAL( KIND=RP ) :: gammaPlus1Div2     , gammaMinus1Div2sg, gammaMinus1Div2g
     REAL( KIND=RP ) :: InvGammaPlus1Div2  , InvGammaMinus1   , InvGamma
     REAL( KIND=RP ) :: gammaDivGammaMinus1, gammaM2 !! = gamma*mach**2
!
!    ========
     CONTAINS
!    ========
!
!     ///////////////////////////////////////////////////////
!
!     --------------------------------------------------
!!    Constructor: Define default values for the physics
!!    variables.
!     --------------------------------------------------
!
      SUBROUTINE ConstructPhysicsStorage( machArg, REArg, PRArg, AOAThetaArg, AOAPhiArg, flowIsNavierStokesArg )
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP) :: machArg, REArg, PRArg
      REAL(KIND=RP) :: AOAPhiArg, AOAThetaArg
      LOGICAL       :: flowIsNavierStokesArg
      
      mach               = machArg
      RE                 = REArg
      PR                 = PRArg
      AOAPhi             = AOAPhiArg
      AOATheta           = AOAThetaArg
      flowIsNavierStokes = flowIsNavierStokesArg
!
      TRef            = 520.0_RP
      TScale          = 198.6_RP
      TRatio          = TScale/TRef
      
      gamma                = 1.4_RP
      gammaMinus1          = gamma - 1.0_RP
      sqrtGamma            = SQRT( gamma )
      gammaMinus1Div2      = gammaMinus1/2.0_RP
      gammaPlus1Div2       = ( gamma + 1.0_RP )/2.0_RP
      gammaMinus1Div2sg    = gammaMinus1Div2 / sqrtGamma
      gammaMinus1Div2g     = gammaMinus1Div2 / gamma
      InvGammaPlus1Div2    = 1.0_RP / gammaPlus1Div2
      InvGammaMinus1       = 1.0_RP / gammaMinus1
      InvGamma             = 1.0_RP / gamma
      gammaDivGammaMinus1  = gamma / gammaMinus1
      gammaM2              = gamma*mach**2
!
      END SUBROUTINE ConstructPhysicsStorage
!
!     ///////////////////////////////////////////////////////
!
!     -------------------------------------------------
!!    Destructor: Does nothing for this storage
!     -------------------------------------------------
!
      SUBROUTINE DestructPhysicsStorage
      
      END SUBROUTINE DestructPhysicsStorage
!
!    **********       
     END MODULE PhysicsStorage
!    **********
!@mark -
!
!  **************
   Module Physics 
!  **************
!
      USE SMConstants
      USE PhysicsStorage
      IMPLICIT NONE
!
!     ---------
!     Constants
!     ---------
!
      INTEGER, PARAMETER   :: WALL_BC = 1, RADIATION_BC = 2
      INTEGER, PARAMETER   :: ROE = 0, LXF = 1, RUSANOV = 2
      REAL(KIND=RP)        :: waveSpeed
      INTEGER              :: boundaryCondition(4), bcType
      INTEGER              :: riemannSolverChoice = RUSANOV
!
!     ========
      CONTAINS 
!     ========
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE RiemannSolver( QLeft, QRight, nHat, flux )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), DIMENSION(N_EQN)  :: Qleft, Qright, flux
         REAL(KIND=RP), DIMENSION(3)      :: nHat
         
         SELECT CASE ( riemannSolverChoice )
            CASE ( ROE )
               CALL RoeSolver( QLeft, QRight, nHat, flux )
            CASE (LXF)
               PRINT *, "3D LXF to be implemented..."
               STOP !DEBUG
               CALL LxFSolver( QLeft, QRight, nHat, flux )
            CASE (RUSANOV)
               CALL RusanovSolver( QLeft, QRight, nHat, flux )               
            CASE DEFAULT
               PRINT *, "Undefined choice of Riemann Solver. Abort"
               STOP
         END SELECT

      
      END SUBROUTINE RiemannSolver
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
         REAL(KIND=RP), DIMENSION(N_EQN) :: Qleft, Qright, flux
         REAL(KIND=RP), DIMENSION(3)     :: nHat
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
      
         rho  = Qleft(1)
         rhou = Qleft(2)
         rhov = Qleft(3)
         rhow = Qleft(4)
         rhoe = Qleft(5)
   
         rhon  = Qright(1)
         rhoun = Qright(2)
         rhovn = Qright(3)
         rhown = Qright(4)
         rhoen = Qright(5)
   
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
         REAL(KIND=RP), DIMENSION(N_EQN) :: Qleft, Qright, flux
         REAL(KIND=RP), DIMENSION(3)     :: nHat
         REAL(KIND=RP)                   :: ds = 1.0_RP
!
!        ---------------
!        Local Variables
!        ---------------
!
!
      REAL(KIND=RP) :: rho , rhou , rhov  , rhoe
      REAL(KIND=RP) :: rhon, rhoun, rhovn , rhoen
      REAL(KIND=RP) :: ul  , vl   , pleft , ql, cl
      REAL(KIND=RP) :: ur  , vr   , pright, qr, cr
      REAL(KIND=RP) :: sM
      REAL(KIND=RP), DIMENSION(N_EQN) :: FL, FR
      
      rho  = Qleft(1)
      rhou = Qleft(2)
      rhov = Qleft(3)
      rhoe = Qleft(4)

      rhon  = Qright(1)
      rhoun = Qright(2)
      rhovn = Qright(3)
      rhoen = Qright(4)

      ul = rhou/rho 
      vl = rhov/rho 
      pleft = (gamma-1.d0)*(rhoe - 0.5d0/rho*(rhou**2 + rhov**2)) 
!
      ur = rhoun/rhon 
      vr = rhovn/rhon 
      pright = (gamma-1.d0)*(rhoen - 0.5d0/rhon*(rhoun**2 + rhovn**2)) 
!
      ql = nHat(1)*ul + nHat(2)*vl 
      qr = nHat(1)*ur + nHat(2)*vr 
      cl = SQRT( gamma*pleft/rho )
      cr = SQRT( gamma*pright/rhon )
!
      FL(1) = rho*ql
      FL(2) = rhou*ql + pleft*nHat(1)
      FL(3) = rhov*ql + pleft*nHat(2)
      FL(4) = (rhoe + pleft)*ql
!
      FR(1) = rhon*qr
      FR(2) = rhoun*qr + pright*nHat(1)
      FR(3) = rhovn*qr + pright*nHat(2)
      FR(4) = (rhoen + pright)*qr
!
      sM = MAX( ABS(ql) + cl, ABS(qr) + cr )
!
      flux = ds * 0.5_RP * ( FL + FR - sM*(Qright - Qleft) )      
         
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
         REAL(KIND=RP), DIMENSION(N_EQN) :: Qleft, Qright, flux
         REAL(KIND=RP), DIMENSION(3)     :: nHat
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
         !REAL(KIND=RP) :: gamma = 1.4_RP
      
         rho  = Qleft(1)
         rhou = Qleft(2)
         rhov = Qleft(3)
         rhow = Qleft(4)
         rhoe = Qleft(5)
   
         rhon  = Qright(1)
         rhoun = Qright(2)
         rhovn = Qright(3)
         rhown = Qright(4)
         rhoen = Qright(5)
   
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
         ar = SQRT(ar2)
         al = SQRT(al2)
!           
         rql = rho*ql 
         rqr = rhon*qr             
         flux(1) = ds*(rql + rqr) 
         flux(2) = ds*(rql*ul + pleft*nHat(1) + rqr*ur + pright*nHat(1)) 
         flux(3) = ds*(rql*vl + pleft*nHat(2) + rqr*vr + pright*nHat(2))
         flux(4) = ds*(rql*wl + pleft*nHat(3) + rqr*wr + pright*nHat(3)) 
         flux(5) = ds*(rql*hl + rqr*hr) 

         smax = MAX(ar+ABS(qr),al+ABS(ql))

         flux = (flux - ds*smax*(Qright-Qleft))/2.d0

         RETURN 
         
      END SUBROUTINE RusanovSolver           
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE xFlux( Q, f )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), DIMENSION(N_EQN) :: Q
         REAL(KIND=RP), DIMENSION(N_EQN) :: f
!
!        ---------------
!        Local Variables
!        ---------------
!
         REAL(KIND=RP) :: u, v, w, rho, rhou, rhov, rhoe, rhow, p
!      
         rho  = Q(1)
         rhou = Q(2)
         rhov = Q(3)
         rhow = Q(4)
         rhoe = Q(5)
!
         u = rhou/rho 
         v = rhov/rho
         w = rhow/rho
         p = gammaMinus1*(rhoe - 0.5_RP*rho*(u**2 + v**2 + w**2)) 
!
         f(1) = rhou 
         f(2) = p + rhou*u 
         f(3) = rhou*v 
         f(4) = rhou*w 
         f(5) = u*(rhoe + p) 
         
      END SUBROUTINE xFlux
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE yFlux( Q, g )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), DIMENSION(N_EQN) :: Q
         REAL(KIND=RP), DIMENSION(N_EQN) :: g
!
!        ---------------
!        Local Variables
!        ---------------
!
         REAL(KIND=RP) :: u, v, w, rho, rhou, rhov, rhoe, rhow, p
!      
         rho  = Q(1)
         rhou = Q(2)
         rhov = Q(3)
         rhow = Q(4)
         rhoe = Q(5)
!
         u = rhou/rho 
         v = rhov/rho 
         w = rhow/rho
         p = gammaMinus1*(rhoe - 0.5_RP*rho*(u**2 + v**2 + w**2)) 
!
         g(1) = rhov 
         g(2) = rhou*v 
         g(3) = p + rhov*v 
         g(4) = rhow*v 
         g(5) = v*(rhoe + p) 
         
      END SUBROUTINE yFlux
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE zFlux( Q, h )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), DIMENSION(N_EQN) :: Q
         REAL(KIND=RP), DIMENSION(N_EQN) :: h
!
!        ---------------
!        Local Variables
!        ---------------
!
         REAL(KIND=RP) :: u, v, w, rho, rhou, rhov, rhoe, rhow, p
!      
         rho  = Q(1)
         rhou = Q(2)
         rhov = Q(3)
         rhow = Q(4)
         rhoe = Q(5)
!
         u = rhou/rho 
         v = rhov/rho 
         w = rhow/rho
         p = gammaMinus1*(rhoe - 0.5_RP*rho*(u**2 + v**2 + w**2)) 
!
         h(1) = rhow 
         h(2) = rhou*w 
         h(3) = rhov*w
         h(4) = p + rhow*w 
         h(5) = w*(rhoe + p) 
         
      END SUBROUTINE zFlux
!
! /////////////////////////////////////////////////////////////////////
!
!@mark -
!---------------------------------------------------------------------
!! DiffusionRiemannSolution computes the coupling on the solution for
!! the calculation of the gradient terms.
!---------------------------------------------------------------------
!
      SUBROUTINE DiffusionRiemannSolution( nHat, QLeft, QRight, Q )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), DIMENSION(N_EQN) :: Qleft, Qright, Q
         REAL(KIND=RP), DIMENSION(3)     :: nHat
!
!        ---------------
!        Local Variables
!        ---------------
!
         INTEGER :: j
!
!        -----------------------------------------------
!        For now, this is simply the Bassi/Rebay average
!        -----------------------------------------------
!
         DO j = 1, N_EQN
            Q(j) = 0.5_RP*(Qleft(j) + Qright(j))
         END DO

      END SUBROUTINE DiffusionRiemannSolution
!
! /////////////////////////////////////////////////////////////////////////////
!
!-----------------------------------------------------------------------------
!! DiffusionRiemannSolution computes the coupling on the gradients for
!! the calculation of the contravariant diffusive flux.
!-----------------------------------------------------------------------------
!
      SUBROUTINE DiffusionRiemannFlux(nHat, ds, Q, gradLeft, gradRight, flux)
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(N_EQN)   :: Q,flux
      REAL(KIND=RP), DIMENSION(3,N_EQN) :: gradLeft, gradRight
      REAL(KIND=RP), DIMENSION(3)       :: nHat
      REAL(KIND=RP)                     :: ds
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER                           :: j,k
      REAL(KIND=RP), DIMENSION(3,N_EQN) :: grad
      REAL(KIND=RP), DIMENSION(N_EQN)   :: fx, fy, fz
!
!     -------------------------------------------------
!     For now, this simply uses the Bassi/Rebay average
!     -------------------------------------------------
!
      DO j = 1, N_GRAD_EQN
         DO k = 1,3
            grad(k,j) = 0.5_RP*(gradLeft(k,j) + gradRight(k,j))
         END DO
      END DO
!
!     ----------------------------
!     Compute the component fluxes
!     ----------------------------
!
      CALL xDiffusiveFlux( Q, grad, fx )
      CALL yDiffusiveFlux( Q, grad, fy )
      CALL zDiffusiveFlux( Q, grad, fz )
!
!     ------------------------------
!     Compute the contravariant flux
!     ------------------------------
!
      DO j = 1, N_EQN
         flux(j) = ds*(nHat(1)*fx(j) + nHat(2)*fy(j) + nHat(3)*fz(j))
      END DO
      
      END SUBROUTINE DiffusionRiemannFlux
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! xDiffusiveFlux computes the x viscous flux component.
!---------------------------------------------------------------------
!
      SUBROUTINE xDiffusiveFlux( Q, grad, f )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
!!    Q contains the solution values
!
      REAL(KIND=RP), DIMENSION(N_EQN)      :: Q
!
!!    grad contains the (physical) gradients needed for the
!!    equations. For the Navier-Stokes equations these are
!!    grad(u), grad(v), grad(w), grad(T).
!
      REAL(KIND=RP), DIMENSION(3,N_GRAD_EQN) :: grad
!
!!     f is the viscous flux in the physical x direction returned by
!!     this routine.
! 
      REAL(KIND=RP), DIMENSION(N_EQN)      :: f
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: tauXX, tauXY, tauXZ
      REAL(KIND=RP) :: T, muOfT, kappaOfT, divVelocity
      REAL(KIND=RP) :: u, v, w
!      
      T        = Temperature(Q)
      muOfT    = MolecularDiffusivity(T)
      kappaOfT = ThermalDiffusivity(T)
      u        = Q(2)/Q(1)
      v        = Q(3)/Q(1)
      w        = Q(4)/Q(1)
      
      divVelocity = grad(1,1) + grad(2,2) + grad(3,3)
      tauXX       = 2.0_RP*muOfT*(grad(1,1) - divVelocity/3._RP)
      tauXY       = muOfT*(grad(1,2) + grad(2,1))
      tauXZ       = muOfT*(grad(1,3) + grad(3,1))
      
      f(1) = 0.0_RP
      f(2) = tauXX/RE
      f(3) = tauXY/RE
      f(4) = tauXZ/RE
      f(5) = (u*tauXX + v*tauXY + w*tauXZ + &
     &        gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*grad(1,4))/RE

      END SUBROUTINE xDiffusiveFlux
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! yDiffusiveFlux computes the y viscous flux component.
!---------------------------------------------------------------------
!
      SUBROUTINE yDiffusiveFlux( Q, grad, f )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
!!    Q contains the solution values
!
      REAL(KIND=RP), DIMENSION(N_EQN)      :: Q
!
!!    grad contains the (physical) gradients needed for the
!!    equations. For the Navier-Stokes equations these are
!!    grad(u), grad(v), grad(w), grad(T).
!
      REAL(KIND=RP), DIMENSION(3,N_GRAD_EQN) :: grad
!
!!     f is the viscous flux in the physical x direction returned by
!!     this routine.
! 
      REAL(KIND=RP), DIMENSION(N_EQN)      :: f
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: tauYX, tauYY, tauYZ
      REAL(KIND=RP) :: T, muOfT, kappaOfT, divVelocity
      REAL(KIND=RP) :: u, v, w
!      
      T        = Temperature(Q)
      muOfT    = MolecularDiffusivity(T)
      kappaOfT = ThermalDiffusivity(T)
      u        = Q(2)/Q(1)
      v        = Q(3)/Q(1)
      w        = Q(4)/Q(1)
      
      divVelocity = grad(1,1) + grad(2,2) + grad(3,3)
      tauYX       = muOfT*(grad(1,2) + grad(2,1))
      tauYY       = 2.0_RP*muOfT*(grad(2,2) - divVelocity/3._RP)
      tauYZ       = muOfT*(grad(2,3) + grad(3,2))
      
      f(1) = 0.0_RP
      f(2) = tauYX/RE
      f(3) = tauYY/RE
      f(4) = tauYZ/RE
      f(5) = (u*tauYX + v*tauYY + w*tauYZ + &
     &        gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*grad(2,4))/RE

      END SUBROUTINE yDiffusiveFlux
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! yDiffusiveFlux computes the y viscous flux component.
!---------------------------------------------------------------------
!
      SUBROUTINE zDiffusiveFlux( Q, grad, f )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
!!    Q contains the solution values
!
      REAL(KIND=RP), DIMENSION(N_EQN)      :: Q
!
!!    grad contains the (physical) gradients needed for the
!!    equations. For the Navier-Stokes equations these are
!!    grad(u), grad(v), grad(w), grad(T).
!
      REAL(KIND=RP), DIMENSION(3,N_GRAD_EQN) :: grad
!
!!     f is the viscous flux in the physical x direction returned by
!!     this routine.
! 
      REAL(KIND=RP), DIMENSION(N_EQN)      :: f
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP)           :: tauZX, tauZY, tauZZ
      REAL(KIND=RP)           :: T, muOfT, kappaOfT, divVelocity
      REAL(KIND=RP)           :: u, v, w
!      
      T        = Temperature(Q)
      muOfT    = MolecularDiffusivity(T)
      kappaOfT = ThermalDiffusivity(T)
      u        = Q(2)/Q(1)
      v        = Q(3)/Q(1)
      w        = Q(4)/Q(1)
      
      divVelocity = grad(1,1) + grad(2,2) + grad(3,3)
      tauZX       = muOfT*(grad(1,3) + grad(3,1))
      tauZY       = muOfT*(grad(2,3) + grad(3,2))
      tauZZ       = 2.0_RP*muOfT*(grad(3,3) - divVelocity/3._RP)
      
      f(1) = 0.0_RP
      f(2) = tauZX/RE
      f(3) = tauZY/RE
      f(4) = tauZZ/RE
      f(5) = (u*tauZX + v*tauZY + w*tauZZ + &
     &        gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*grad(3,4))/RE

      END SUBROUTINE zDiffusiveFlux
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! GradientValuesForQ takes the solution (Q) values and returns the
!! quantities of which the gradients will be taken.
!---------------------------------------------------------------------
!
      SUBROUTINE GradientValuesForQ( Q, U )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(N_EQN)     , INTENT(IN)  :: Q
      REAL(KIND=RP), DIMENSION(N_GRAD_EQN), INTENT(OUT) :: U
!
!     ---------------
!     Local Variables
!     ---------------
!     
      U(1) = Q(2)/Q(1)
      U(2) = Q(3)/Q(1)
      U(3) = Q(4)/Q(1)
      U(4) = Temperature(Q)

      END SUBROUTINE GradientValuesForQ
!
! /////////////////////////////////////////////////////////////////////
!
!@mark -
!---------------------------------------------------------------------
!! Compute the pressure from the state variables
!---------------------------------------------------------------------
!
      FUNCTION Pressure(Q) RESULT(P)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(N_EQN) :: Q
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: P
      
      P = gammaMinus1*(Q(5) - 0.5_RP*(Q(2)**2 + Q(3)**2 + Q(4)**2)/Q(1))

      END FUNCTION Pressure
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! Compute the molecular diffusivity by way of Sutherland's law
!---------------------------------------------------------------------
!
      FUNCTION MolecularDiffusivity(T) RESULT(mu)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP) :: T !! The temperature
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: mu !! The diffusivity
!      
      mu = (1._RP + tRatio)/(T + tRatio)*T*SQRT(T)


      END FUNCTION MolecularDiffusivity
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! Compute the thermal diffusivity by way of Sutherland's law
!---------------------------------------------------------------------
!
      FUNCTION ThermalDiffusivity(T) RESULT(kappa)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP) :: T !! The temperature
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: kappa !! The diffusivity
!      
      kappa = (1._RP + tRatio)/(T + tRatio)*T*SQRT(T)


      END FUNCTION ThermalDiffusivity
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! Compute the temperature from the state variables
!---------------------------------------------------------------------
!
      FUNCTION Temperature(Q) RESULT(T)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(N_EQN) :: Q
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: T
!
      T = gammaM2*Pressure(Q)/Q(1)

      END FUNCTION Temperature
      
   END Module Physics
!
! /////////////////////////////////////////////////////////////////////
!
!----------------------------------------------------------------------
!! This routine returns the maximum eigenvalues for the Euler equations 
!! for the given solution value in each spatial direction. 
!! These are to be used to compute the local time step.
!----------------------------------------------------------------------
!
      SUBROUTINE ComputeEigenvaluesForState( Q, eigen )
      
      USE SMConstants
      USE PhysicsStorage
      USE Physics, ONLY:Pressure
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=Rp), DIMENSION(N_EQN) :: Q
      REAL(KIND=Rp), DIMENSION(3)     :: eigen
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=Rp) :: u, v, w, p, a
!      
      u = ABS( Q(2)/Q(1) )
      v = ABS( Q(3)/Q(1) )
      w = ABS( Q(4)/Q(1) )
      p = Pressure(Q)
      a = SQRT(gamma*p/Q(1))
      
      eigen(1) = u + a
      eigen(2) = v + a
      eigen(3) = w + a
      
      END SUBROUTINE ComputeEigenvaluesForState
