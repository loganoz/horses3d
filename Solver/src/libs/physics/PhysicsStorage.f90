!
!//////////////////////////////////////////////////////
!
!   @File:    PhysicsStorage.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Wed Dec  6 17:42:24 2017
!   @Last revision date: Tue Dec 26 20:56:54 2017
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 11dd5b683eba468bc0e8e4db146cbc1fbc952a8a
!
!//////////////////////////////////////////////////////
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
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: LAMBDA_STABILIZATION_KEY  = "lambda stabilization"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: COMPUTE_GRADIENTS_KEY     = "compute gradients"
         
         CHARACTER(LEN=KEYWORD_LENGTH), DIMENSION(2) :: physicsKeywords = [MACH_NUMBER_KEY, FLOW_EQUATIONS_KEY]
         
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: CENTRAL_SOLVER_NAME      ="central"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: ROE_SOLVER_NAME          ="roe"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: RUSANOV_SOLVER_NAME      ="rusanov"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: LAXFRIEDRICHS_SOLVER_NAME="lax-friedrichs"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: STDROE_SOLVER_NAME       ="standard roe"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: ROEPIKE_SOLVER_NAME      ="roe-pike"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: LOWDISSROE_SOLVER_NAME   ="low dissipation roe"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: MATRIXDISS_SOLVER_NAME   ="matrix dissipation"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: VISCOUSNS_SOLVER_NAME    ="viscous ns"
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

     private
     public    flowIsNavierStokes, N_EQN, N_GRAD_EQN, NDIM, IX, IY, IZ
     public    NCONS, IRHO, IRHOU, IRHOV, IRHOW, IRHOE
     public    NGRAD, IGU, IGV, IGW, IGT
     public    NPRIM, IPIRHO, IPU, IPV, IPW, IPP, IPT, IPA2
     public    TScale, TRatio
     public    Thermodynamics, RefValues, Dimensionless
     public    Thermodynamics_t, RefValues_t, Dimensionless_t
     public    lambdaStab, computeGradients, whichRiemannSolver, whichAverage
     public    RIEMANN_ROE, RIEMANN_LXF, RIEMANN_RUSANOV, RIEMANN_STDROE
     public    RIEMANN_CENTRAL, RIEMANN_ROEPIKE, RIEMANN_LOWDISSROE
     public    RIEMANN_VISCOUSNS, RIEMANN_MATRIXDISS
     public    STANDARD_SPLIT, DUCROS_SPLIT, MORINISHI_SPLIT
     public    KENNEDYGRUBER_SPLIT, PIROZZOLI_SPLIT, ENTROPYCONS_SPLIT
     public    ENTROPYANDENERGYCONS_SPLIT

     public    ConstructPhysicsStorage, DestructPhysicsStorage, DescribePhysicsStorage
     public    CheckPhysicsInputIntegrity
!
!    ----------------------------
!    Either NavierStokes or Euler
!    ----------------------------
!
     logical, protected :: flowIsNavierStokes = .true.
     logical, protected :: computeGradients   = .true.
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
!    ----------------------------------------
!!   The positions of the primitive variables
!    ----------------------------------------
!
     INTEGER, PARAMETER       :: NPRIM = 7
     INTEGER, PARAMETER       :: IPIRHO = 1, IPU = 2, IPV = 3, IPW = 4, IPP = 5, IPT = 6, IPA2 = 7
!
!    ---------------------------------------
!!   The positions of the gradient variables
!    ---------------------------------------
!
     INTEGER, PARAMETER  :: NGRAD = 4
     INTEGER, PARAMETER  :: IGU = 1 , IGV = 2 , IGW = 3 , IGT = 4
!
!    --------------------------------------------
!    The temperature scale in the Sutherland law:
!    198.6 for temperatures in R, 110.3 for
!    temperatures in K.
!    --------------------------------------------
!
     REAL( KIND=RP ), protected :: TScale
!
!    ------------------------------------------------
!    The ratio of the scale and reference tempartures
!    ------------------------------------------------
!
     REAL( KIND=RP ), protected :: TRatio 
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
     integer            :: whichAverage               = -1
!
!    -------------------------------------
!    Lambda stabilization - 1.0 by default
!    -------------------------------------
!
     real(kind=RP), protected :: lambdaStab = 1.0_RP     
!
!    ------------------------------
!    Thermodynamic specs of the Air
!    ------------------------------
!
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
              287.15_RP * 5.0_RP / 9.0_RP / ( 1.4_RP - 1.0_RP ), & ! cv
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
      SUBROUTINE ConstructPhysicsStorage( controlVariables, success )
      USE FTValueDictionaryClass
      USE PhysicsKeywordsModule
      use Utilities, only: toLower
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(FTValueDictionary) :: controlVariables
      LOGICAL                 :: success
!
!     ---------------
!     Local variables
!     ---------------
!
      CHARACTER(LEN=KEYWORD_LENGTH) :: keyword
      type(Thermodynamics_t), pointer  :: thermodynamics_
      type(RefValues_t)                :: refValues_
      type(Dimensionless_t)            :: dimensionless_
!
!     --------------------
!     Collect input values
!     --------------------
!
      success = .TRUE.
      CALL CheckPhysicsInputIntegrity(controlVariables,success)
      IF(.NOT. success) RETURN 
!
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
      dimensionless_ % Mach = controlVariables % doublePrecisionValueForKey(MACH_NUMBER_KEY)
      dimensionless_ % Pr   = 0.72_RP
      dimensionless_ % Fr   = 0.0_RP
!
!     *********************
!     Select flow equations
!     *********************
!
      keyword = controlVariables % stringValueForKey(FLOW_EQUATIONS_KEY,KEYWORD_LENGTH)
      CALL toLower(keyword)

      IF ( keyword == "euler" )     THEN
         flowIsNavierStokes = .FALSE.
         dimensionless_ % Re = 0.0_RP   
         dimensionless_ % mu = 0.0_RP
         dimensionless_ % kappa = 0.0_RP

      ELSE 
         flowIsNavierStokes = .TRUE.
!
!        ----------------------------
!        Look for the Reynolds number
!        ----------------------------
!
         IF ( controlVariables % containsKey(REYNOLDS_NUMBER_KEY) )     THEN
            dimensionless_ % Re = controlVariables % doublePrecisionValueForKey(REYNOLDS_NUMBER_KEY) 

         ELSE 
            PRINT *, "Input file is missing entry for keyword: ", REYNOLDS_NUMBER_KEY
            success = .FALSE.
            RETURN 

         END IF 
!
!        ------------------------------------------------
!        Set molecular viscosity and thermal conductivity
!        ------------------------------------------------
!
         dimensionless_ % mu   = 1.0_RP / dimensionless_ % Re
         dimensionless_ % kappa = 1.0_RP / ( thermodynamics_ % gammaMinus1 * &
                                              POW2( dimensionless_ % Mach) * &
                                      dimensionless_ % Re * dimensionless_ % Pr )
      END IF 
!
!     **************************************
!     Check if state gradients are requested
!     **************************************
!
      if ( controlVariables % containsKey(COMPUTE_GRADIENTS_KEY) ) then
!
!        Do not compute gradients if Euler equations are used and it is specified in the control file
!        --------------------------------------------------------------------------------------------
         if ( .not. flowIsNavierStokes ) then
            if ( .not. controlVariables % logicalValueForKey(COMPUTE_GRADIENTS_KEY) ) then
               computeGradients = .false.

            end if
         end if
      else
!
!        Do not compute gradients if Euler equations and the default option is selected
!        ------------------------------------------------------------------------------
         if ( .not. flowIsNavierStokes ) then
            computeGradients = .false.

         end if
      end if

      dimensionless_ % gammaM2 = thermodynamics_ % gamma * POW2( dimensionless_ % Mach )
      dimensionless_ % invFroudeSquare = 0.0_RP
!
!     ********************
!     Set reference values: TODO read from parameter file
!                           Ok, but be sure to change the mesh reading accordingly (x = x / refValues % L)
!     ********************
!
      refValues_ % L = 1.0_RP       ! m
      refValues_ % T = 520.0_RP     ! Rankine

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
!
!     *********************************************
!     Choose the Riemann solver (by default is Roe)
!     *********************************************
!
      IF ( controlVariables % containsKey(RIEMANN_SOLVER_NAME_KEY) )     THEN
!
!        Get keyword from control variables
!        ----------------------------------
         keyword = controlVariables % stringValueForKey(key             = RIEMANN_SOLVER_NAME_KEY,&
                                                        requestedLength = KEYWORD_LENGTH)
         CALL toLower(keyword)
!
!        Choose the appropriate Riemann Solver
!        -------------------------------------
         select case ( keyword )
         case(ROE_SOLVER_NAME) 
            whichRiemannSolver = RIEMANN_ROE

         case(LAXFRIEDRICHS_SOLVER_NAME)
            whichRiemannSolver = RIEMANN_LXF 

         case(RUSANOV_SOLVER_NAME)
            whichRiemannSolver = RIEMANN_RUSANOV
   
         case(STDROE_SOLVER_NAME)
            whichRiemannSolver = RIEMANN_STDROE

         case(CENTRAL_SOLVER_NAME)
            whichRiemannSolver = RIEMANN_CENTRAL

         case(ROEPIKE_SOLVER_NAME)
            whichRiemannSolver = RIEMANN_ROEPIKE

         case(LOWDISSROE_SOLVER_NAME)
            whichRiemannSolver = RIEMANN_LOWDISSROE

         case(MATRIXDISS_SOLVER_NAME)
            whichRiemannSolver = RIEMANN_MATRIXDISS

         case(VISCOUSNS_SOLVER_NAME)
            whichRiemannSolver = RIEMANN_VISCOUSNS
            
            
         case default 
            print*, "Riemann solver: ", trim(keyword), " is not implemented."
            print*, "Options available are:"
            print*, "   * Central"
            print*, "   * Roe"
            print*, "   * Standard Roe"
            print*, "   * Lax-Friedrichs"
            print*, "   * Rusanov"
            print*, "   * Roe-Pike"
            print*, "   * Low dissipation Roe"
            print*, "   * Matrix dissipation"
            print*, "   * Viscous NS"
            errorMessage(STD_OUT)
            stop
         end select 
      else 
!
!        Select Roe by default
!        ---------------------
         whichRiemannSolver = RIEMANN_ROE 

      end if
!
!     --------------------
!     Lambda stabilization
!     --------------------
!
      if ( controlVariables % containsKey(LAMBDA_STABILIZATION_KEY)) then
         lambdaStab = controlVariables % doublePrecisionValueForKey(LAMBDA_STABILIZATION_KEY)

      else
!
!        By default, lambda is 1 (full upwind stabilization)
!        ---------------------------------------------------
         lambdaStab = 1.0_RP

      end if
!
!     --------------------------------------------------
!     If central fluxes are used, set lambdaStab to zero
!     --------------------------------------------------
!
      if ( whichRiemannSolver .eq. RIEMANN_CENTRAL ) lambdaStab = 0.0_RP
!
!     ***************
!     Angle of attack
!     ***************
!
      IF ( controlVariables % containsKey(AOA_PHI_KEY) )     THEN
         refValues_ % AOAPhi = controlVariables % doublePrecisionValueForKey(AOA_PHI_KEY) 

      ELSE
!
!        Phi angle of attack is zero by default
!        --------------------------------------
         refValues_ % AOAPhi = 0.0_RP

      END IF 

      IF ( controlVariables % containsKey(AOA_THETA_KEY) )     THEN
         refValues_ % AOATheta = controlVariables % doublePrecisionValueForKey(AOA_THETA_KEY) 

      ELSE
!
!        Theta angle of attack is zero by default
!        ----------------------------------------
         refValues_ % AOATheta = 0.0_RP

      END IF 
!
!     **************************
!     Sutherland's law constants
!     **************************
!
      TScale          = 198.6_RP
      TRatio          = TScale/ refValues_ % T
!
!     **********************************************************************
!     Set the global (proteted) thermodynamics, dimensionless, and refValues
!     **********************************************************************
!
      call setThermodynamics( thermodynamics_ )
      call setDimensionless( dimensionless_ )
      call setRefValues( refValues_ )
!
!     ********
!     Describe
!     ********
!
      CALL DescribePhysicsStorage()

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
         use MPI_Process_Info
         IMPLICIT NONE
         real(kind=RP)  :: pRef

         if ( .not. MPI_Process % isRoot ) return 

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

