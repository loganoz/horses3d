!      Physics.f90
!      Created: 2011-07-20 09:17:26 -0400 
!      By: David Kopriva
!      From DSEM Code
!
!    Special version for Test problem with F = Q, averaging for riemann
!    solver.
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
     REAL( KIND=RP ) :: AOA
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
      SUBROUTINE ConstructPhysicsStorage( machArg, REArg, PRArg, flowIsNavierStokesArg )
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP) :: machArg, REArg, PRArg
      LOGICAL       :: flowIsNavierStokesArg
      
      mach               = machArg
      RE                 = REArg
      PR                 = PRArg
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
      INTEGER, PARAMETER   :: ROE = 0, LXF = 1
      REAL(KIND=RP)        :: waveSpeed
      INTEGER              :: boundaryCondition(4), bcType
      INTEGER              :: riemannSolverChoice = ROE
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
!
!        ----------------------------------------------------------------
!        F = Q\hat x + Q\hat y + Q \hat z$ so the normal flux is as below
!        ----------------------------------------------------------------
!
         flux = 0.5_RP*(Qleft + Qright)*( nHat(1) + nHat(2) + nHat(3))
      
      END SUBROUTINE RiemannSolver
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

         f = Q
                  
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

         g = Q
                  
         
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

         h = Q
         
      END SUBROUTINE zFlux
!
! /////////////////////////////////////////////////////////////////////
!
!@mark -
!---------------------------------------------------------------------
!! DiffusionRiemannSolution comutes the coupling on the solution for
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
!! DiffusionRiemannSolution comutes the coupling on the gradients for
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
      p = gammaMinus1*(Q(4) - 0.5_Rp*(Q(2)**2 + Q(3)**2 + Q(4)**2 )/Q(1))
      a = SQRT(gamma*p/Q(1))
      
      eigen(1) = u + a
      eigen(2) = v + a
      eigen(3) = w + a
      
      END SUBROUTINE ComputeEigenvaluesForState
