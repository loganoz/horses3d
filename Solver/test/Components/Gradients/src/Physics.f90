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
#include "Includes.h"
      Module PhysicsKeywordsModule
         IMPLICIT NONE 
         INTEGER, PARAMETER :: KEYWORD_LENGTH = 132
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: MACH_NUMBER_KEY           = "mach number"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: REYNOLDS_NUMBER_KEY       = "reynolds number"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PRANDTL_NUMBER_KEY        = "prandtl number"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: FROUDE_NUMBER_KEY         = "froude number"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: AOA_THETA_KEY             = "aoa theta"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: AOA_PHI_KEY               = "aoa phi"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: FLOW_EQUATIONS_KEY        = "flow equations"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: RIEMANN_SOLVER_NAME_KEY   = "riemann solver"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: LESMODEL_KEY              = "les model"
         
         CHARACTER(LEN=KEYWORD_LENGTH), DIMENSION(2) :: physicsKeywords = [MACH_NUMBER_KEY, FLOW_EQUATIONS_KEY]
         
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: ROE_SOLVER_NAME           = "roe"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: RUSANOV_SOLVER_NAME       = "rusanov"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: LAXFRIEDRICHS_SOLVER_NAME = "lax friedrichs"

         !PARTICLES 
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: particlesKey             = "lagrangian particles"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: numberOfParticlesKey     = "number of particles"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: STOKES_NUMBER_PART_KEY   = "stokes number"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: GAMMA_PART_KEY           = "gamma"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PHI_M_PART_KEY           = "phi_m"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: I0_PART_KEY              = "radiation source"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: gx_PART_KEY              = "gravity_x"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: gy_PART_KEY              = "gravity_y"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: gz_PART_KEY              = "gravity_z"

         
      END MODULE PhysicsKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
!    
!    ******
     MODULE PhysicsStorage
!    ******
!
     USE SMConstants
     use FluidData
     
     IMPLICIT NONE
     SAVE
!
!    ----------------------------
!    Either NavierStokes or Euler
!    ----------------------------
!
     LOGICAL :: flowIsNavierStokes = .true.
     LOGICAL :: computeGradients = .true.
     logical :: useLESModel = .false.
!
!    --------------------------
!!   The sizes of the NS system
!    --------------------------
!
     INTEGER, PARAMETER :: N_EQN = 5, N_GRAD_EQN = 4
!
!    -------------------------------------------
!!   The positions of the conservative variables
!    -------------------------------------------
!
     INTEGER, PARAMETER       :: NCONS = 5
     INTEGER, PARAMETER       :: IRHO = 1 , IRHOU = 2 , IRHOV = 3 , IRHOW = 4 , IRHOE = 5
!
!    ---------------------------------------
!!   The positions of the gradient variables
!    ---------------------------------------
!
     INTEGER, PARAMETER  :: IGU = 1 , IGV = 2 , IGW = 3 , IGT = 4
     integer, parameter  :: NGRAD = 4
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
     real(kind=RP)   :: Lref, timeRef
!    ----------------------------------
!
!    --------------------------
!    Riemann solver definitions
!    --------------------------
!
     integer, parameter :: RIEMANN_ROE        = 0
     integer, parameter :: RIEMANN_LXF        = 1
     integer, parameter :: RIEMANN_RUSANOV    = 2
     integer, parameter :: RIEMANN_STDROE     = 4
     integer, parameter :: RIEMANN_CENTRAL    = 5
     integer, parameter :: RIEMANN_ROEPIKE    = 6
     integer, parameter :: RIEMANN_LOWDISSROE = 7
     integer, parameter :: RIEMANN_VISCOUSNS  = 8
     integer, parameter :: RIEMANN_MATRIXDISS = 9
     integer, protected :: whichRiemannSolver = -1
     real(kind=RP)      :: lambdaStab = 0.0_RP
!
!    -----------------------------
!    Available averaging functions
!    -----------------------------
!
     integer, parameter :: STANDARD_SPLIT             = 1
     integer, parameter :: MORINISHI_SPLIT            = 2
     integer, parameter :: DUCROS_SPLIT               = 3
     integer, parameter :: KENNEDYGRUBER_SPLIT        = 4
     integer, parameter :: PIROZZOLI_SPLIT            = 5
     integer, parameter :: ENTROPYCONS_SPLIT          = 6
     integer, parameter :: ENTROPYANDENERGYCONS_SPLIT = 7
     integer            :: whichAverage = -1


     type(Thermodynamics_t), target, private :: ThermodynamicsAir = Thermodynamics_t( &
                                                              "Air", & ! Name
                                    287.15_RP * 5.0_RP / 9.0_RP, & ! R
                                                         1.4_RP, & ! gamma
                                                   sqrt(1.4_RP), & ! sqrtGamma
                                                1.4_RP - 1.0_RP, & ! gammaMinus1         
                                     (1.4_RP - 1.0_RP) / 2.0_RP, & ! gammaMinus1Div2
                                     (1.4_RP + 1.0_RP) / 2.0_RP, & ! gammaPlus1Div2
                    (1.4_RP - 1.0_RP) / (2.0_RP * sqrt(1.4_RP)), & ! gammaMinus1Div2sg 
                          (1.4_RP - 1.0_RP) / (2.0_RP * 1.4_RP), & ! gammaMinus1Div2g 
                                     2.0_RP / (1.4_RP + 1.0_RP), & ! InvGammaPlus1Div2 
                                     1.0_RP / (1.4_RP - 1.0_RP), & ! InvGammaMinus1
                                                1.0_RP / 1.4_RP, & ! InvGamma
                                   1.4_RP / ( 1.4_RP - 1.0_RP ), & ! gammaDivGammaMinus1
     287.15_RP * 5.0_RP / 9.0_RP * 1.4_RP / ( 1.4_RP - 1.0_RP ), & ! cp
              287.15_RP * 5.0_RP / 9.0_RP / ( 1.4_RP - 1.0_RP ), & ! cp
                                                         0.0_RP  & ! Bulk viscosity ratio
)
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
!
!     ---------------
!     Local variables
!     ---------------
!
      type(Thermodynamics_t), pointer  :: thermodynamics_
      type(RefValues_t)                :: refValues_
      type(Dimensionless_t)            :: dimensionless_
!
!     ---------------------
!     Set the gas to be air
!     ---------------------
!
      thermodynamics_ => thermodynamicsAir
!
!     ------------------------
!     Dimensionless quantities
!     ------------------------
!
      dimensionless_ % cp = thermodynamics_ % gamma * thermodynamics_ % InvGammaMinus1
      dimensionless_ % cv = thermodynamics_ % InvGammaMinus1
      dimensionless_ % Mach = machArg
      dimensionless_ % Pr   = PRArg
      dimensionless_ % Re  = ReArg
      flowIsNavierStokes = flowIsNavierStokesArg

      if ( flowIsNavierStokes ) then
         dimensionless_ % mu   = 1.0_RP / dimensionless_ % Re
         dimensionless_ % kappa = 1.0_RP / ( thermodynamics_ % gammaMinus1 * &
                                              POW2( dimensionless_ % Mach) * &
                                      dimensionless_ % Re * dimensionless_ % Pr )
      else
         dimensionless_ % mu = 0.0_RP
         dimensionless_ % kappa = 0.0_RP
      end if
      dimensionless_ % gammaM2 = thermodynamics_ % gamma * POW2( dimensionless_ % Mach )

!      if ( dimensionless_ % Fr == huge(1.d0) ) then  
!            dimensionless_ % invFroudeSquare = 0.0_RP 
!      else  
!            dimensionless_ % invFroudeSquare = 1.0_RP / POW2( dimensionless_ % Fr ) 
!      endif  
!
!     ----------------
!     Reference values
!     ----------------
!
      Lref  = 1.0_RP 
      refValues_ % T  = 520.0_RP
      refValues_ % rho = 101325.0_RP / (thermodynamics_ % R * refValues_ % T)
      refValues_ % V =   dimensionless_ % Mach &
                       * sqrt( thermodynamics_ % gamma * thermodynamics_ % R * refValues_ % T )
      refValues_ % p = refValues_ % rho * POW2( refValues_ % V )

      if ( flowIsNavierStokes ) then
         refValues_ % mu = refValues_ % rho * refValues_ % V * Lref / dimensionless_ % Re
         refValues_ % kappa = refValues_ % mu * thermodynamics_ % cp / dimensionless_ % Pr

      else
         refValues_ % mu = 0.0_RP
         refValues_ % kappa = 0.0_RP
      
      end if

      timeref = Lref / refValues_ % V

      TScale          = 198.6_RP
      TRatio          = TScale/ refValues_ % T

      call setThermodynamics(thermodynamics_)
      call setDimensionless(dimensionless_)
      call setRefValues(refValues_)
      
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
!     //////////////////////////////////////////////////////
!
!     -----------------------------------------
!!    Descriptor: Shows the gathered data
!     -----------------------------------------
!
      SUBROUTINE DescribePhysicsStorage()
         USE Headers
         IMPLICIT NONE
         real(kind=RP)  :: pRef

         pRef = thermodynamics % R * refValues % rho * refValues % T

         write(STD_OUT,'(/,/)')
         if (flowIsNavierStokes) then
            call Section_Header("Loading Navier-Stokes physics")
         else
            call Section_Header("Loading Euler physics")
         end if

         write(STD_OUT,'(/)')
         call SubSection_Header("Fluid data")
         write(STD_OUT,'(30X,A,A22,A10)') "->" , "Gas: " , "Air"
         write(STD_OUT,'(30X,A,A22,F10.3,A)') "->" , "State constant: " , thermodynamics % R, " I.S."
         write(STD_OUT,'(30X,A,A22,F10.3)') "->" , "Specific heat ratio: " , thermodynamics % gamma

         write(STD_OUT,'(/)')
         call SubSection_Header("Reference quantities")
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference Temperature: " , refValues % T, " K."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference pressure: " , pRef, " Pa."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference density: " , refValues % rho , " kg/m^3."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference velocity: " , refValues % V , " m/s."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reynolds length: " , Lref , " m."
         
         if ( flowIsNavierStokes ) then
            write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference viscosity: ",refValues % mu , " Pa·s."
            write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference conductivity: ", refValues % kappa, " W/(m·K)."
         end if

         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference time: ", timeref, " s."

         write(STD_OUT,'(/)')
         call SubSection_Header("Dimensionless quantities")
         write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Mach number: " , dimensionless % Mach
         if ( flowIsNavierStokes ) then
            write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Reynolds number: " , dimensionless % Re
            write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Prandtl number: " , dimensionless % Pr
         end if

         write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Froude number: " , dimensionless % Fr 

      END SUBROUTINE DescribePhysicsStorage
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE CheckPhysicsInputIntegrity( controlVariables, success )  
         USE FTValueDictionaryClass
         USE PhysicsKeywordsModule
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(FTValueDictionary) :: controlVariables
         LOGICAL                 :: success
!
!        ---------------
!        Local variables
!        ---------------
!
         CLASS(FTObject), POINTER :: obj
         INTEGER                  :: i
         success = .TRUE.
         
         DO i = 1, SIZE(physicsKeywords)
            obj => controlVariables % objectForKey(physicsKeywords(i))
            IF ( .NOT. ASSOCIATED(obj) )     THEN
               PRINT *, "Input file is missing entry for keyword: ",physicsKeywords(i)
               success = .FALSE. 
            END IF  
         END DO  
         
      END SUBROUTINE CheckPhysicsInputIntegrity
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
      REAL(KIND=RP)        :: waveSpeed
      INTEGER              :: boundaryCondition(4), bcType


!
!    ---------------
!    Interface block
!    ---------------
!
     interface InviscidFlux
         module procedure InviscidFlux0D , InviscidFlux3D
     end interface InviscidFlux

     interface ViscousFlux
         module procedure ViscousFlux0D , ViscousFlux2D, ViscousFlux3D
         module procedure ViscousFlux0DWithSGS , ViscousFlux3DWithSGS
     end interface ViscousFlux
    
     interface SutherlandsLaw
         module procedure MolecularDiffusivity
     end interface SutherlandsLaw
!
!     ========
      CONTAINS 
!     ========
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      pure subroutine InviscidFlux0D( Q, F )
         implicit none
         real(kind=RP), intent(in)           :: Q(1:NCONS)
         real(kind=RP), intent(out)          :: F(1:NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)           :: u , v , w , p

         F(:,IX) = Q
         F(:,IY) = Q
         F(:,IZ) = Q

      end subroutine InviscidFlux0D

      pure subroutine InviscidFlux3D( N , Q, F )
         implicit none
         integer,       intent (in) :: N(3)
         real(kind=RP), intent (in) :: Q(1:NCONS, 0:N(1) , 0:N(2) , 0:N(3))
         real(kind=RP), intent(out) :: F(1:NCONS, 0:N(1) , 0:N(2) , 0:N(3), 1:NDIM)

         F(:,:,:,:,IX) = Q
         F(:,:,:,:,IY) = Q
         F(:,:,:,:,IZ) = Q

      end subroutine InviscidFlux3D
!
!     -------------------------------------------------------------------------------
!     Subroutine for computing the Jacobian of the inviscid flux when it has the form 
!
!        F = f*iHat + g*jHat + h*kHat
!
!     First index indicates the flux term and second index indicates the conserved 
!     variable term. 
!     ***** This routine is necessary for computing the analytical Jacobian. *****
!     -------------------------------------------------------------------------------
      pure subroutine InviscidJacobian(q,dfdq,dgdq,dhdq)
         implicit none
         !-------------------------------------------------
         real(kind=RP), intent (in)  :: q(NCONS)
         real(kind=RP), intent (out) :: dfdq(NCONS,NCONS)
         real(kind=RP), intent (out) :: dgdq(NCONS,NCONS)
         real(kind=RP), intent (out) :: dhdq(NCONS,NCONS)
         !-------------------------------------------------
         real(kind=RP)  :: u,v,w ! Velocity components
         real(kind=RP)  :: V2    ! Total velocity squared
         real(kind=RP)  :: p     ! Pressure
         real(kind=RP)  :: H     ! Total enthalpy
         !-------------------------------------------------
         
         associate( gammaMinus1 => thermodynamics % gammaMinus1, & 
                    gamma => thermodynamics % gamma )
         
         u  = q(IRHOU) / q(IRHO)
         v  = q(IRHOV) / q(IRHO)
         w  = q(IRHOW) / q(IRHO)
         V2 = u*u + v*v + w*w
!         p  = Pressure(q)
         H  = (q(IRHOE) + p) / q(IRHO)
!
!        Flux in the x direction (f)
!        ---------------------------

         dfdq(1,1) = 0._RP
         dfdq(1,2) = 1._RP
         dfdq(1,3) = 0._RP
         dfdq(1,4) = 0._RP
         dfdq(1,5) = 0._RP
         
         dfdq(2,1) = -u*u + 0.5_RP*gammaMinus1*V2
         dfdq(2,2) = (3._RP - gamma) * u
         dfdq(2,3) = -gammaMinus1 * v
         dfdq(2,4) = -gammaMinus1 * w
         dfdq(2,5) = gammaMinus1
         
         dfdq(3,1) = -u*v
         dfdq(3,2) = v
         dfdq(3,3) = u
         dfdq(3,4) = 0._RP
         dfdq(3,5) = 0._RP
         
         dfdq(4,1) = -u*w
         dfdq(4,2) = w
         dfdq(4,3) = 0._RP
         dfdq(4,4) = u
         dfdq(4,5) = 0._RP
         
         dfdq(5,1) = u * (0.5_RP*gammaMinus1*V2 - H)
         dfdq(5,2) = H - gammaMinus1 * u*u
         dfdq(5,3) = -gammaMinus1 * u*v
         dfdq(5,4) = -gammaMinus1 * u*w
         dfdq(5,5) = gamma * u
         
!
!        Flux in the y direction (g)
!        ---------------------------
         
         dgdq(1,1) = 0._RP
         dgdq(1,2) = 0._RP
         dgdq(1,3) = 1._RP
         dgdq(1,4) = 0._RP
         dgdq(1,5) = 0._RP
         
         dgdq(2,1) = -u*v
         dgdq(2,2) = v
         dgdq(2,3) = u
         dgdq(2,4) = 0._RP
         dgdq(2,5) = 0._RP
         
         dgdq(3,1) = -v*v + 0.5_RP*gammaMinus1*V2
         dgdq(3,2) = -gammaMinus1 * u
         dgdq(3,3) = (3._RP - gamma) * v
         dgdq(3,4) = -gammaMinus1 * w
         dgdq(3,5) = gammaMinus1
         
         dgdq(4,1) = -v*w
         dgdq(4,2) = 0._RP
         dgdq(4,3) = w
         dgdq(4,4) = v
         dgdq(4,5) = 0._RP
         
         dgdq(5,1) = v * (0.5_RP*gammaMinus1*V2 - H)
         dgdq(5,2) = -gammaMinus1 * u*v
         dgdq(5,3) = H - gammaMinus1 * v*v
         dgdq(5,4) = -gammaMinus1 * v*w
         dgdq(5,5) = gamma * v
!
!        Flux in the z direction (h)
!        ---------------------------
         
         dhdq(1,1) = 0._RP
         dhdq(1,2) = 0._RP
         dhdq(1,3) = 0._RP
         dhdq(1,4) = 1._RP
         dhdq(1,5) = 0._RP
         
         dhdq(2,1) = -u*w
         dhdq(2,2) = w
         dhdq(2,3) = 0._RP
         dhdq(2,4) = u
         dhdq(2,5) = 0._RP
         
         dhdq(3,1) = -v*w
         dhdq(3,2) = 0._RP
         dhdq(3,3) = w
         dhdq(3,4) = v
         dhdq(3,5) = 0._RP
         
         dhdq(4,1) = -w*w + 0.5_RP*gammaMinus1*V2
         dhdq(4,2) = -gammaMinus1 * u
         dhdq(4,3) = -gammaMinus1 * v
         dhdq(4,4) = (3._RP - gamma) * w
         dhdq(4,5) = gammaMinus1
         
         dhdq(5,1) = w * (0.5_RP*gammaMinus1*V2 - H)
         dhdq(5,2) = -gammaMinus1 * u*w
         dhdq(5,3) = -gammaMinus1 * v*w
         dhdq(5,4) = H - gammaMinus1 * w*w
         dhdq(5,5) = gamma * w
         
         end associate
         
      end subroutine InviscidJacobian
!
! /////////////////////////////////////////////////////////////////////
!
!@mark -
!---------------------------------------------------------------------
!! DiffusionRiemannSolution computes the coupling on the solution for
!! the calculation of the gradient terms.
!---------------------------------------------------------------------
!
      pure subroutine ViscousFlux0DWithSGS( Q , U_x , U_y , U_z, mu, kappa, tauSGS, qSGS, F)
         implicit none
         real ( kind=RP ) , intent ( in ) :: Q    ( 1:NCONS          ) 
         real ( kind=RP ) , intent ( in ) :: U_x  ( 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_y  ( 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_z  ( 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: mu, kappa
         real(kind=RP),    intent(in)     :: tauSGS(NDIM,NDIM), qSGS(NDIM)
         real(kind=RP), intent(out)       :: F    ( 1:NCONS , 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                    :: T , muOfT , kappaOfT
         real(kind=RP)                    :: divV
         real(kind=RP)                    :: u , v , w

         F = 0.0_RP
         F(1:N_GRAD_EQN,IX) = U_x
         F(1:N_GRAD_EQN,IY) = U_y
         F(1:N_GRAD_EQN,IZ) = U_z

      end subroutine ViscousFlux0DWithSGS

      pure subroutine ViscousFlux3DWithSGS( N, Q , U_x , U_y , U_z, mu, kappa, tauSGS, qSGS, F)
         implicit none
         integer          , intent ( in ) :: N(3)
         real ( kind=RP ) , intent ( in ) :: Q    ( 1:NCONS, 0:N(1) , 0:N(2) , 0:N(3)) 
         real ( kind=RP ) , intent ( in ) :: U_x  ( 1:N_GRAD_EQN, 0:N(1) , 0:N(2) , 0:N(3)) 
         real ( kind=RP ) , intent ( in ) :: U_y  ( 1:N_GRAD_EQN, 0:N(1) , 0:N(2) , 0:N(3)) 
         real ( kind=RP ) , intent ( in ) :: U_z  ( 1:N_GRAD_EQN, 0:N(1) , 0:N(2) , 0:N(3)) 
         real ( kind=RP ) , intent ( in ) :: mu   ( 0:N(1) , 0:N(2) , 0:N(3)) 
         real ( kind=RP ) , intent ( in ) :: kappa   ( 0:N(1) , 0:N(2) , 0:N(3)) 
         real ( kind=RP ) , intent ( in ) :: tauSGS  (NDIM, NDIM, 0:N(1) , 0:N(2) , 0:N(3)) 
         real ( kind=RP ) , intent ( in ) :: qSGS    (NDIM, 0:N(1) , 0:N(2) , 0:N(3)) 
         real ( kind=RP ) , intent ( out) :: F    ( 1:NCONS, 0:N(1) , 0:N(2) , 0:N(3), 1:NDIM )

         F = 0.0_RP
         F(1:N_GRAD_EQN,:,:,:,IX) = U_x
         F(1:N_GRAD_EQN,:,:,:,IY) = U_y
         F(1:N_GRAD_EQN,:,:,:,IZ) = U_z

      end subroutine ViscousFlux3DWithSGS

      pure subroutine ViscousFlux0D( Q , U_x , U_y , U_z, mu, kappa, F)
         implicit none
         real ( kind=RP ) , intent ( in ) :: Q    ( 1:NCONS          ) 
         real ( kind=RP ) , intent ( in ) :: U_x  ( 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_y  ( 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_z  ( 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: mu, kappa
         real(kind=RP), intent(out)       :: F    ( 1:NCONS , 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                    :: T , muOfT , kappaOfT
         real(kind=RP)                    :: divV
         real(kind=RP)                    :: u , v , w

         F = 0.0_RP
         F(1:N_GRAD_EQN,IX) = U_x
         F(1:N_GRAD_EQN,IY) = U_y
         F(1:N_GRAD_EQN,IZ) = U_z

      end subroutine ViscousFlux0D

      pure subroutine ViscousFlux2D( N, Q , U_x , U_y , U_z, mu, kappa, F)
         implicit none
         integer          , intent ( in ) :: N(2)
         real ( kind=RP ) , intent ( in ) :: Q    ( 1:NCONS, 0:N(1) , 0:N(2)) 
         real ( kind=RP ) , intent ( in ) :: U_x  ( 1:N_GRAD_EQN, 0:N(1) , 0:N(2)) 
         real ( kind=RP ) , intent ( in ) :: U_y  ( 1:N_GRAD_EQN, 0:N(1) , 0:N(2)) 
         real ( kind=RP ) , intent ( in ) :: U_z  ( 1:N_GRAD_EQN, 0:N(1) , 0:N(2)) 
         real ( kind=RP ) , intent ( in ) :: mu   ( 0:N(1) , 0:N(2)) 
         real ( kind=RP ) , intent ( in ) :: kappa   ( 0:N(1) , 0:N(2)) 
         real ( kind=RP ) , intent ( out) :: F    ( 1:NCONS, 1:NDIM, 0:N(1) , 0:N(2))

         F = 0.0_RP
         F(1:N_GRAD_EQN,IX,:,:) = U_x
         F(1:N_GRAD_EQN,IY,:,:) = U_y
         F(1:N_GRAD_EQN,IZ,:,:) = U_z

      end subroutine ViscousFlux2D

      pure subroutine ViscousFlux3D( N, Q , U_x , U_y , U_z, mu, kappa, F)
         implicit none
         integer          , intent ( in ) :: N(3)
         real ( kind=RP ) , intent ( in ) :: Q    ( 1:NCONS, 0:N(1) , 0:N(2) , 0:N(3)) 
         real ( kind=RP ) , intent ( in ) :: U_x  ( 1:N_GRAD_EQN, 0:N(1) , 0:N(2) , 0:N(3)) 
         real ( kind=RP ) , intent ( in ) :: U_y  ( 1:N_GRAD_EQN, 0:N(1) , 0:N(2) , 0:N(3)) 
         real ( kind=RP ) , intent ( in ) :: U_z  ( 1:N_GRAD_EQN, 0:N(1) , 0:N(2) , 0:N(3)) 
         real ( kind=RP ) , intent ( in ) :: mu   ( 0:N(1) , 0:N(2) , 0:N(3)) 
         real ( kind=RP ) , intent ( in ) :: kappa   ( 0:N(1) , 0:N(2) , 0:N(3)) 
         real ( kind=RP ) , intent ( out) :: F    ( 1:NCONS, 0:N(1) , 0:N(2) , 0:N(3), 1:NDIM )

         F = 0.0_RP
         F(1:N_GRAD_EQN,:,:,:,IX) = U_x
         F(1:N_GRAD_EQN,:,:,:,IY) = U_y
         F(1:N_GRAD_EQN,:,:,:,IZ) = U_z

      end subroutine ViscousFlux3D
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! Compute the molecular diffusivity by way of Sutherland's law
!---------------------------------------------------------------------
!
      PURE FUNCTION MolecularDiffusivity(T) RESULT(mu)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), INTENT(IN) :: T !! The temperature
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
      PURE FUNCTION ThermalDiffusivity(T) RESULT(kappa)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), INTENT(IN) :: T !! The temperature
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: kappa !! The diffusivity
!      
      kappa = (1._RP + tRatio)/(T + tRatio)*T*SQRT(T)


      END FUNCTION ThermalDiffusivity

     pure subroutine getStressTensor(Q,U_x,U_y,U_z,tau)
         implicit none
         real(kind=RP), intent(in)      :: Q   (1:NCONS         )
         real(kind=RP), intent(in)      :: U_x (1:N_GRAD_EQN    )
         real(kind=RP), intent(in)      :: U_y (1:N_GRAD_EQN    )
         real(kind=RP), intent(in)      :: U_z (1:N_GRAD_EQN    )
         real(kind=RP), intent(out)     :: tau (1:NDIM, 1:NDIM   )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: T , muOfT
         real(kind=RP) :: divV

         associate ( mu0 => dimensionless % mu )

!         T     = Temperature(Q)
!         muOfT = SutherlandsLaw(T)

         divV = U_x(IGU) + U_y(IGV) + U_z(IGW)

         tau(IX,IX) = mu0 * muOfT * (2.0_RP * U_x(IGU) - 2.0_RP/3.0_RP * divV )
         tau(IY,IX) = mu0 * muOfT * ( U_x(IGV) + U_y(IGU) )
         tau(IZ,IX) = mu0 * muOfT * ( U_x(IGW) + U_z(IGU) )
         tau(IX,IY) = tau(IY,IX)
         tau(IY,IY) = mu0 * muOfT * (2.0_RP * U_y(IGV) - 2.0_RP/3.0_RP * divV )
         tau(IZ,IY) = mu0 * muOfT * ( U_y(IGW) + U_z(IGV) )
         tau(IX,IZ) = tau(IZ,IX)
         tau(IY,IZ) = tau(IZ,IY)
         tau(IZ,IZ) = mu0 * muOfT * (2.0_RP * U_z(IGW) - 2.0_RP/3.0_RP * divV )

         end associate

      end subroutine getStressTensor


!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! Compute the temperature from the state variables
!---------------------------------------------------------------------
!

      
   END Module Physics
   
   module RiemannSolvers
      use SMConstants
      use Physics
      use PhysicsStorage
      use FluidData, only: equationOfState, getThermalConductivity
      contains
         subroutine SetRiemannSolver(which, splitType)
            integer, intent(in)  :: which
            integer, intent(in)  :: splitType
         end subroutine SetRiemannSolver
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE RiemannSolver( QLeft, QRight, nHat, t1, t2, flux )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), DIMENSION(N_EQN)  :: Qleft, Qright, flux
         REAL(KIND=RP), DIMENSION(3)      :: nHat, t1, t2
         
         flux = 0.5_RP*(Qleft + Qright)*( nHat(1) + nHat(2) + nHat(3) )
      
      END SUBROUTINE RiemannSolver
      
      SUBROUTINE RiemannSolver_dFdQ(ql,qr,nHat,dfdq_num,side)
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), DIMENSION(N_EQN)        :: Ql, Qr
         REAL(KIND=RP), DIMENSION(N_EQN,N_EQN)  :: dfdq_num
         integer                                :: side
         REAL(KIND=RP), DIMENSION(3)            :: nHat
         
      
      END SUBROUTINE RiemannSolver_dFdQ
   end module RiemannSolvers
!@mark -
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
!      USE Physics, ONLY:Pressure
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
!      p = Pressure(Q)
      a = SQRT(thermodynamics % gamma*p/Q(1))
      
      eigen(1) = u + a
      eigen(2) = v + a
      eigen(3) = w + a
      
      END SUBROUTINE ComputeEigenvaluesForState
