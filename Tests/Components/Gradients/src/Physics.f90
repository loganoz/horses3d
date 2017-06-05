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
     REAL( KIND=RP ) :: PR = 0.72_RP
!
!    ----------------------------------------
!!   The free-stream or reference temperature
!!   with default in R.
!    ----------------------------------------
!
     REAL( KIND=RP ) :: TRef 
!
!    ----------------------------------------
!!   The free-stream or reference pressure
!!   with default in Pa.
!    ----------------------------------------
!
     REAL( KIND=RP ) :: pRef 
!
!    ----------------------------------------
!!   The length in the Reynolds number
!    ----------------------------------------
!
     REAL( KIND=RP ) :: reynoldsLength 
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
!    ----------------------
!!   The gas state constant
!    ----------------------
!
     REAL( KIND=RP ) :: Rgas
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
!    ------------------------------------
!    Riemann solver associated quantities
!    ------------------------------------
!
     INTEGER, PARAMETER :: ROE = 0, LXF = 1, RUSANOV = 2
     INTEGER            :: riemannSolverChoice = ROE
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
!     //////////////////////////////////////////////////////
!
!     -----------------------------------------
!!    Descriptor: Shows the gathered data
!     -----------------------------------------
!
      SUBROUTINE DescribePhysicsStorage()
         USE Headers
         IMPLICIT NONE

         write(STD_OUT,'(/,/)')
         if (flowIsNavierStokes) then
            call Section_Header("Loading Navier-Stokes physics")
         else
            call Section_Header("Loading Euler physics")
         end if

         write(STD_OUT,'(/)')
         call SubSection_Header("Fluid data")
         write(STD_OUT,'(30X,A,A22,A10)') "->" , "Gas: " , "Air"
         write(STD_OUT,'(30X,A,A22,F10.3,A)') "->" , "State constant: " , Rgas, " I.S."
         write(STD_OUT,'(30X,A,A22,F10.3)') "->" , "Specific heat ratio: " , gamma

         write(STD_OUT,'(/)')
         call SubSection_Header("Reference quantities")
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference Temperature: " , TRef, " K."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference pressure: " , pRef, " Pa."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference density: " , pRef / (Rgas * TRef) , " kg/m^3."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference velocity: " , Mach * sqrt(gamma * Rgas * TRef) , " m/s."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reynolds length: " , reynoldsLength , " m."
         
         if ( flowIsNavierStokes ) then
            write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference viscosity: ", &
                     sqrt(gamma) * Mach * reynoldsLength * pRef / ( RE * sqrt(Rgas * TRef) ) , " Pa·s."
            write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference conductivity: ", &
                     gammaDivGammaMinus1 * Rgas * sqrt(gamma) * Mach * reynoldsLength * pRef / ( RE * sqrt(Rgas * TRef) ) / PR, &
                     " W/(m·K)."
         end if
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference time: " , &
                     reynoldsLength / (Mach * sqrt(gamma * Rgas * TRef) ) , " s."

         write(STD_OUT,'(/)')
         call SubSection_Header("Dimensionless quantities")
         write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Mach number: " , Mach
         if ( flowIsNavierStokes ) then
            write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Reynolds number: " , RE
            write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Prandtl number: " , PR
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
         real(kind=RP), intent (in) :: Q(0:N , 1:NCONS)
         real(kind=RP)              :: F(0:N , 1:NCONS , 1:NDIM)
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
         real(kind=RP), intent (in) :: Q(0:N , 0:N , 1:NCONS)
         real(kind=RP)              :: F(0:N , 0:N , 1:NCONS , 1:NDIM)
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

      pure function InviscidFlux3D( N , Q ) result ( F )
         implicit none
         integer,       intent (in) :: N
         real(kind=RP), intent (in) :: Q(0:N , 0:N , 0:N , 1:NCONS)
         real(kind=RP)              :: F(0:N , 0:N , 0:N , 1:NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)           :: u(0:N,0:N,0:N) , v(0:N,0:N,0:N) , w(0:N,0:N,0:N) , p(0:N,0:N,0:N)

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
         real ( kind=RP ) , intent ( in ) :: Q    ( 0:N , 1:NCONS          ) 
         real ( kind=RP ) , intent ( in ) :: U_x  ( 0:N , 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_y  ( 0:N , 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_z  ( 0:N , 1:N_GRAD_EQN     ) 
         real(kind=RP)                    :: F    ( 0:N , 1:NCONS , 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                    :: T(0:N) , muOfT(0:N) , kappaOfT(0:N)
         real(kind=RP)                    :: divV(0:N)
         real(kind=RP)                    :: u(0:N) , v(0:N) , w(0:N)
         integer                          :: i


         F = 0.0_RP
         F(:,1:N_GRAD_EQN,IX) = U_x
         F(:,1:N_GRAD_EQN,IY) = U_y
         F(:,1:N_GRAD_EQN,IZ) = U_z

      end function ViscousFlux1D

      pure function ViscousFlux2D( N , Q , U_x , U_y , U_z ) result (F)
         implicit none
         integer          , intent ( in ) :: N
         real ( kind=RP ) , intent ( in ) :: Q    ( 0:N , 0:N , 1:NCONS          ) 
         real ( kind=RP ) , intent ( in ) :: U_x  ( 0:N , 0:N , 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_y  ( 0:N , 0:N , 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_z  ( 0:N , 0:N , 1:N_GRAD_EQN     ) 
         real(kind=RP)                    :: F    ( 0:N , 0:N , 1:NCONS , 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                    :: T(0:N,0:N) , muOfT(0:N,0:N) , kappaOfT(0:N,0:N)
         real(kind=RP)                    :: divV(0:N,0:N)
         real(kind=RP)                    :: u(0:N,0:N) , v(0:N,0:N) , w(0:N,0:N)
         integer                          :: i , j 


         F = 0.0_RP
         F(:,:,1:N_GRAD_EQN,IX) = U_x
         F(:,:,1:N_GRAD_EQN,IY) = U_y
         F(:,:,1:N_GRAD_EQN,IZ) = U_z


      end function ViscousFlux2D

      pure function ViscousFlux3D( N , Q , U_x , U_y , U_z ) result (F)
         implicit none
         integer          , intent ( in ) :: N
         real ( kind=RP ) , intent ( in ) :: Q    ( 0:N , 0:N , 0:N , 1:NCONS          ) 
         real ( kind=RP ) , intent ( in ) :: U_x  ( 0:N , 0:N , 0:N , 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_y  ( 0:N , 0:N , 0:N , 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_z  ( 0:N , 0:N , 0:N , 1:N_GRAD_EQN     ) 
         real ( kind=RP )                 :: F    ( 0:N , 0:N , 0:N , 1:NCONS , 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: T(0:N,0:N,0:N) , muOfT(0:N,0:N,0:N) , kappaOfT(0:N,0:N,0:N)
         real(kind=RP) :: divV(0:N,0:N,0:N)
         real(kind=RP) :: u(0:N,0:N,0:N) , v(0:N,0:N,0:N) , w(0:N,0:N,0:N)
         integer       :: i , j , k


         F = 0.0_RP
         F(:,:,:,1:N_GRAD_EQN,IX) = U_x
         F(:,:,:,1:N_GRAD_EQN,IY) = U_y
         F(:,:,:,1:N_GRAD_EQN,IZ) = U_z


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

      SUBROUTINE GradientValuesForQ_3D( Q, U )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), INTENT(IN)  :: Q(0:,0:,0:,:)
      REAL(KIND=RP), INTENT(OUT) :: U(0:,0:,0:,:)
      integer                    :: N 
!
!     ---------------
!     Local Variables
!     ---------------
!     
      U(:,:,:,1) = Q(:,:,:,1)
      U(:,:,:,2) = Q(:,:,:,2)
      U(:,:,:,3) = Q(:,:,:,3)
      U(:,:,:,4) = Q(:,:,:,4)

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
      
      P = gammaMinus1*(Q(5) - 0.5_RP*(Q(2)**2 + Q(3)**2 + Q(4)**2)/Q(1))

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
      T = gammaM2*Pressure(Q)/Q(1)

      END FUNCTION Temperature
      
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
      a = SQRT(gamma*p/Q(1))
      
      eigen(1) = u + a
      eigen(2) = v + a
      eigen(3) = w + a
      
      END SUBROUTINE ComputeEigenvaluesForState
