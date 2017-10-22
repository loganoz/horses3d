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
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: AOA_THETA_KEY             = "aoa theta"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: AOA_PHI_KEY               = "aoa phi"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: FLOW_EQUATIONS_KEY        = "flow equations"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: RIEMANN_SOLVER_NAME_KEY   = "riemann solver"
         
         CHARACTER(LEN=KEYWORD_LENGTH), DIMENSION(2) :: physicsKeywords = [MACH_NUMBER_KEY, FLOW_EQUATIONS_KEY]
         
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: ROE_SOLVER_NAME           = "roe"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: RUSANOV_SOLVER_NAME       = "rusanov"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: LAXFRIEDRICHS_SOLVER_NAME = "lax friedrichs"
         
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
!
!    --------------------------
!!   The sizes of the NS system
!    --------------------------
!
     INTEGER, PARAMETER :: N_EQN = 5, N_GRAD_EQN = 4
!
!    -----------------------------
!    Number of physical dimensions
!    -----------------------------
!
     INTEGER, PARAMETER       :: NDIM = 3
     INTEGER, PARAMETER       :: IX = 1 , IY = 2 , IZ = 3
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
!    ----------------------------------
!
!    ------------------------------------
!    Riemann solver associated quantities
!    ------------------------------------
!
     INTEGER, PARAMETER :: ROE = 0, LXF = 1, RUSANOV = 2
     INTEGER            :: riemannSolverChoice = ROE

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
      dimensionless_ % invFroudeSquare = 0.0_RP
!
!     ----------------
!     Reference values
!     ----------------
!
      refValues_ % L  = 1.0_RP 
      refValues_ % T  = 520.0_RP
      refValues_ % rho = 101325.0_RP / (thermodynamics_ % R * refValues_ % T)
      refValues_ % V =   dimensionless_ % Mach &
                       * sqrt( thermodynamics_ % gamma * thermodynamics_ % R * refValues_ % T )
      refValues_ % p = refValues_ % rho * POW2( refValues_ % V )

      if ( flowIsNavierStokes ) then
         refValues_ % mu = refValues_ % rho * refValues_ % V * refValues_ % L / dimensionless_ % Re
         refValues_ % kappa = refValues_ % mu * thermodynamics_ % cp / dimensionless_ % Pr

      else
         refValues_ % mu = 0.0_RP
         refValues_ % kappa = 0.0_RP
      
      end if

      refValues_ % time = refValues_ % L / refValues_ % V

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
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reynolds length: " , refValues % L , " m."
         
         if ( flowIsNavierStokes ) then
            write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference viscosity: ",refValues % mu , " Pa·s."
            write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference conductivity: ", refValues % kappa, " W/(m·K)."
         end if

         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference time: ", refValues % time, " s."

         write(STD_OUT,'(/)')
         call SubSection_Header("Dimensionless quantities")
         write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Mach number: " , dimensionless % Mach
         if ( flowIsNavierStokes ) then
            write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Reynolds number: " , dimensionless % Re
            write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Prandtl number: " , dimensionless % Pr
         end if

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
     interface GradientValuesForQ
         module procedure GradientValuesForQ_0D , GradientValuesForQ_3D
     end interface GradientValuesForQ

     interface InviscidFlux
         module procedure InviscidFlux0D , InviscidFlux1D , InviscidFlux2D , InviscidFlux3D
     end interface InviscidFlux

     interface ViscousFlux
         module procedure ViscousFlux0D , ViscousFlux1D , ViscousFlux2D , ViscousFlux3D
     end interface ViscousFlux
    
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
         
         flux = 0.5_RP*(Qleft + Qright)*( nHat(1) + nHat(2) + nHat(3) )
      
      END SUBROUTINE RiemannSolver
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      pure function InviscidFlux0D( Q ) result ( F )
         implicit none
         real(kind=RP), intent(in)           :: Q(1:NCONS)
         real(kind=RP)           :: F(1:NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)           :: u , v , w , p

         F(:,IX) = Q
         F(:,IY) = Q
         F(:,IZ) = Q

      end function InviscidFlux0D

      pure function InviscidFlux1D( N , Q ) result ( F )
         implicit none
         integer,       intent (in) :: N
         real(kind=RP), intent (in) :: Q(1:NCONS, 0:N)
         real(kind=RP)              :: F(1:NCONS, 0:N, 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)           :: u(0:N) , v(0:N) , w(0:N) , p(0:N)

         F(:,:,IX) = Q
         F(:,:,IY) = Q
         F(:,:,IZ) = Q 

      end function InviscidFlux1D

      pure function InviscidFlux2D( N , Q ) result ( F )
         implicit none
         integer,       intent (in) :: N
         real(kind=RP), intent (in) :: Q(1:NCONS, 0:N , 0:N)
         real(kind=RP)              :: F(1:NCONS, 0:N , 0:N, 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)           :: u(0:N,0:N) , v(0:N,0:N) , w(0:N,0:N) , p(0:N,0:N)

         F(:,:,:,IX) = Q 
         F(:,:,:,IY) = Q 
         F(:,:,:,IZ) = Q 

      end function InviscidFlux2D

      pure function InviscidFlux3D( Nx, Ny, Nz , Q ) result ( F )
         implicit none
         integer,       intent (in) :: Nx, Ny, Nz
         real(kind=RP), intent (in) :: Q(1:NCONS, 0:Nx , 0:Ny , 0:Nz)
         real(kind=RP)              :: F(1:NCONS, 0:Nx , 0:Ny , 0:Nz, 1:NDIM)

         F(:,:,:,:,IX) = Q
         F(:,:,:,:,IY) = Q
         F(:,:,:,:,IZ) = Q

      end function InviscidFlux3D
!
! /////////////////////////////////////////////////////////////////////
!
!@mark -
!---------------------------------------------------------------------
!! DiffusionRiemannSolution computes the coupling on the solution for
!! the calculation of the gradient terms.
!---------------------------------------------------------------------
!
      pure function ViscousFlux0D( Q , U_x , U_y , U_z ) result (F)
         implicit none
         real ( kind=RP ) , intent ( in ) :: Q    ( 1:NCONS          ) 
         real ( kind=RP ) , intent ( in ) :: U_x  ( 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_y  ( 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_z  ( 1:N_GRAD_EQN     ) 
         real(kind=RP)                    :: F    ( 1:NCONS , 1:NDIM )
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

      end function ViscousFlux0D

      pure function ViscousFlux1D( N , Q , U_x , U_y , U_z ) result (F)
         implicit none
         integer          , intent ( in ) :: N
         real ( kind=RP ) , intent ( in ) :: Q    ( 1:NCONS, 0:N ) 
         real ( kind=RP ) , intent ( in ) :: U_x  ( 1:N_GRAD_EQN, 0:N ) 
         real ( kind=RP ) , intent ( in ) :: U_y  ( 1:N_GRAD_EQN, 0:N ) 
         real ( kind=RP ) , intent ( in ) :: U_z  ( 1:N_GRAD_EQN, 0:N ) 
         real(kind=RP)                    :: F    ( 1:NCONS, 0:N, 1:NDIM )

         F = 0.0_RP
         F(1:N_GRAD_EQN,:,IX) = U_x
         F(1:N_GRAD_EQN,:,IY) = U_y
         F(1:N_GRAD_EQN,:,IZ) = U_z

      end function ViscousFlux1D

      pure function ViscousFlux2D( N , Q , U_x , U_y , U_z ) result (F)
         implicit none
         integer          , intent ( in ) :: N
         real ( kind=RP ) , intent ( in ) :: Q    ( 1:NCONS, 0:N , 0:N ) 
         real ( kind=RP ) , intent ( in ) :: U_x  ( 1:N_GRAD_EQN, 0:N , 0:N ) 
         real ( kind=RP ) , intent ( in ) :: U_y  ( 1:N_GRAD_EQN, 0:N , 0:N ) 
         real ( kind=RP ) , intent ( in ) :: U_z  ( 1:N_GRAD_EQN, 0:N , 0:N ) 
         real(kind=RP)                    :: F    ( 1:NCONS, 0:N , 0:N, 1:NDIM )

         F = 0.0_RP
         F(1:N_GRAD_EQN,:,:,IX) = U_x
         F(1:N_GRAD_EQN,:,:,IY) = U_y
         F(1:N_GRAD_EQN,:,:,IZ) = U_z

      end function ViscousFlux2D

      pure function ViscousFlux3D( Nx, Ny, Nz , Q , U_x , U_y , U_z ) result (F)
         implicit none
         integer          , intent ( in ) :: Nx, Ny, Nz
         real ( kind=RP ) , intent ( in ) :: Q    ( 1:NCONS, 0:Nx , 0:Ny , 0:Nz) 
         real ( kind=RP ) , intent ( in ) :: U_x  ( 1:N_GRAD_EQN, 0:Nx , 0:Ny , 0:Nz) 
         real ( kind=RP ) , intent ( in ) :: U_y  ( 1:N_GRAD_EQN, 0:Nx , 0:Ny , 0:Nz) 
         real ( kind=RP ) , intent ( in ) :: U_z  ( 1:N_GRAD_EQN, 0:Nx , 0:Ny , 0:Nz) 
         real ( kind=RP )                 :: F    ( 1:NCONS, 0:Nx , 0:Ny , 0:Nz, 1:NDIM )

         F = 0.0_RP
         F(1:N_GRAD_EQN,:,:,:,IX) = U_x
         F(1:N_GRAD_EQN,:,:,:,IY) = U_y
         F(1:N_GRAD_EQN,:,:,:,IZ) = U_z

      end function ViscousFlux3D
!
!
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! GradientValuesForQ takes the solution (Q) values and returns the
!! quantities of which the gradients will be taken.
!---------------------------------------------------------------------
!
      SUBROUTINE GradientValuesForQ_0D( Q, U )
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
      U(1) = Q(1)
      U(2) = Q(2)
      U(3) = Q(3)
      U(4) = Q(4)

      END SUBROUTINE GradientValuesForQ_0D

      SUBROUTINE GradientValuesForQ_3D( Nx, Ny, Nz, Q, U )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      integer,       intent(in)  :: Nx, Ny, Nz
      REAL(KIND=RP), INTENT(IN)  :: Q(1:NCONS,0:Nx,0:Ny,0:Nz)
      REAL(KIND=RP), INTENT(OUT) :: U(1:N_GRAD_EQN,0:Nx,0:Ny,0:Nz)
!
!     ---------------
!     Local Variables
!     ---------------
!     
      U(1,:,:,:) = Q(1,:,:,:)
      U(2,:,:,:) = Q(2,:,:,:)
      U(3,:,:,:) = Q(3,:,:,:)
      U(4,:,:,:) = Q(4,:,:,:)

      END SUBROUTINE GradientValuesForQ_3D
!
! /////////////////////////////////////////////////////////////////////
!
!@mark -
!---------------------------------------------------------------------
!! Compute the pressure from the state variables
!---------------------------------------------------------------------
!
      PURE FUNCTION Pressure(Q) RESULT(P)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(N_EQN), INTENT(IN) :: Q
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: P
      
      P = thermodynamics % gammaMinus1*(Q(5) - 0.5_RP*(Q(2)**2 + Q(3)**2 + Q(4)**2)/Q(1))

      END FUNCTION Pressure
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
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! Compute the temperature from the state variables
!---------------------------------------------------------------------
!
      PURE FUNCTION Temperature(Q) RESULT(T)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(N_EQN), INTENT(IN) :: Q
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: T
!
      T = dimensionless % gammaM2*Pressure(Q)/Q(1)

      END FUNCTION Temperature

      function getStressTensor(Q,U_x,U_y,U_z) result(tau)
         implicit none
         real ( kind=RP ) , intent ( in ) :: Q    ( 1:NCONS          ) 
         real ( kind=RP ) , intent ( in ) :: U_x  ( 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_y  ( 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_z  ( 1:N_GRAD_EQN     ) 
         real(kind=RP)                    :: tau  ( 1:NDIM, 1:NDIM   )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: T , muOfT
         real(kind=RP) :: divV

         associate ( mu0 => dimensionless % mu )

         T     = Temperature(Q)
         muOfT = MolecularDiffusivity(T)

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

      end function getStressTensor

      
   END Module Physics
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
      a = SQRT(thermodynamics % gamma*p/Q(1))
      
      eigen(1) = u + a
      eigen(2) = v + a
      eigen(3) = w + a
      
      END SUBROUTINE ComputeEigenvaluesForState
