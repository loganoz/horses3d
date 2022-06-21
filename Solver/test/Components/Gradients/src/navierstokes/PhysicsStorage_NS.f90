!
!//////////////////////////////////////////////////////
!
!   @File:    PhysicsStorage_NS.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Mon May 14 19:03:29 2018
!   @Last revision date: Mon Jul  2 14:17:30 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 7af1f42fb2bc9ea3a0103412145f2a925b4fac5e
!
!//////////////////////////////////////////////////////
!
      
#include "Includes.h"
Module Physics_NSKeywordsModule
         IMPLICIT NONE 
         INTEGER, PARAMETER :: KEYWORD_LENGTH = 132
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: MACH_NUMBER_KEY           = "mach number"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: REYNOLDS_NUMBER_KEY       = "reynolds number"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PRANDTL_NUMBER_KEY        = "prandtl number"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: AOA_THETA_KEY             = "aoa theta"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: AOA_PHI_KEY               = "aoa phi"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: FLOW_EQUATIONS_KEY        = "flow equations"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: RIEMANN_SOLVER_NAME_KEY   = "riemann solver"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: LAMBDA_STABILIZATION_KEY  = "lambda stabilization"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: LESMODEL_KEY              = "les model"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: COMPUTE_GRADIENTS_KEY     = "compute gradients"
         
         CHARACTER(LEN=KEYWORD_LENGTH), DIMENSION(2) :: physics_NSKeywords = [MACH_NUMBER_KEY, FLOW_EQUATIONS_KEY]
         
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: CENTRAL_SOLVER_NAME      ="central"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: ROE_SOLVER_NAME          ="roe"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: RUSANOV_SOLVER_NAME      ="rusanov"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: LAXFRIEDRICHS_SOLVER_NAME="lax-friedrichs"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: STDROE_SOLVER_NAME       ="standard roe"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: ROEPIKE_SOLVER_NAME      ="roe-pike"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: LOWDISSROE_SOLVER_NAME   ="low dissipation roe"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: MATRIXDISS_SOLVER_NAME   ="matrix dissipation"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: VISCOUSNS_SOLVER_NAME    ="viscous ns"

         !PARTICLES 
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: particlesKey             = "lagrangian particles"         
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: numberOfParticlesKey     = "number of particles"          
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: STOKES_NUMBER_PART_KEY   = "stokes number" 
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: FROUDE_NUMBER_KEY        = "froude number"          
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: GAMMA_PART_KEY           = "gamma" 
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PHI_M_PART_KEY           = "phi_m" 
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: I0_PART_KEY              = "radiation source" 
      END MODULE Physics_NSKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
!    
!    ******
     MODULE PhysicsStorage_NS
!    ******
!
     USE SMConstants
     use FluidData_NS
     
     IMPLICIT NONE

     private
     public    flowIsNavierStokes, NCONS, NGRAD
     public    IRHO, IRHOU, IRHOV, IRHOW, IRHOE
     public    IGU, IGV, IGW, IGT
     public    NPRIM, IPIRHO, IPU, IPV, IPW, IPP, IPT, IPA2
     public    TScale, TRatio
     public    lambdaStab, computeGradients, whichRiemannSolver, whichAverage
     public    RIEMANN_ROE, RIEMANN_LXF, RIEMANN_RUSANOV, RIEMANN_STDROE
     public    RIEMANN_CENTRAL, RIEMANN_ROEPIKE, RIEMANN_LOWDISSROE
     public    RIEMANN_VISCOUSNS, RIEMANN_MATRIXDISS
     public    STANDARD_SPLIT, DUCROS_SPLIT, MORINISHI_SPLIT
     public    KENNEDYGRUBER_SPLIT, PIROZZOLI_SPLIT, ENTROPYCONS_SPLIT
     public    ENTROPYANDENERGYCONS_SPLIT
      
     public    ConstructPhysicsStorage_NS, DestructPhysicsStorage_NS, DescribePhysicsStorage_NS
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
     INTEGER, PARAMETER :: NCONS = 5, NGRAD = 4
!
!    -------------------------------------------
!!   The positions of the conservative variables
!    -------------------------------------------
!
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
      SUBROUTINE ConstructPhysicsStorage_NS( machArg, REArg, PRArg, flowIsNavierStokesArg )
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
!            dimensionless_ % invFr2 = 0.0_RP 
!      else  
!            dimensionless_ % invFr2 = 1.0_RP / POW2( dimensionless_ % Fr ) 
!      endif  
!
!     ----------------
!     Reference values
!     ----------------
!
      refValues_ % T  = 520.0_RP
      refValues_ % rho = 101325.0_RP / (thermodynamics_ % R * refValues_ % T)
      refValues_ % V =   dimensionless_ % Mach &
                       * sqrt( thermodynamics_ % gamma * thermodynamics_ % R * refValues_ % T )
      refValues_ % p = refValues_ % rho * POW2( refValues_ % V )

      if ( flowIsNavierStokes ) then
         refValues_ % mu = refValues_ % rho * refValues_ % V * 1.0_RP / dimensionless_ % Re
         refValues_ % kappa = refValues_ % mu * thermodynamics_ % cp / dimensionless_ % Pr

      else
         refValues_ % mu = 0.0_RP
         refValues_ % kappa = 0.0_RP
      
      end if


      TScale          = 198.6_RP
      TRatio          = TScale/ refValues_ % T

      call setThermodynamics(thermodynamics_)
      call setDimensionless(dimensionless_)
      call setRefValues(refValues_)
      
!
      END SUBROUTINE ConstructPhysicsStorage_NS
!
!     ///////////////////////////////////////////////////////
!
!     -------------------------------------------------
!!    Destructor: Does nothing for this storage
!     -------------------------------------------------
!
      SUBROUTINE DestructPhysicsStorage_NS
      
      END SUBROUTINE DestructPhysicsStorage_NS
!
!     //////////////////////////////////////////////////////
!
!     -----------------------------------------
!!    Descriptor: Shows the gathered data
!     -----------------------------------------
!
      SUBROUTINE DescribePhysicsStorage_NS()
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
         
         if ( flowIsNavierStokes ) then
            write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference viscosity: ",refValues % mu , " Pa·s."
            write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference conductivity: ", refValues % kappa, " W/(m·K)."
         end if


         write(STD_OUT,'(/)')
         call SubSection_Header("Dimensionless quantities")
         write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Mach number: " , dimensionless % Mach
         if ( flowIsNavierStokes ) then
            write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Reynolds number: " , dimensionless % Re
            write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Prandtl number: " , dimensionless % Pr
         end if

         write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Froude number: " , dimensionless % Fr 

      END SUBROUTINE DescribePhysicsStorage_NS
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE CheckPhysicsInputIntegrity( controlVariables, success )  
         USE FTValueDictionaryClass
         USE Physics_NSKeywordsModule
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
         
         DO i = 1, SIZE(physics_NSKeywords)
            obj => controlVariables % objectForKey(physics_NSKeywords(i))
            IF ( .NOT. ASSOCIATED(obj) )     THEN
               PRINT *, "Input file is missing entry for keyword: ",physics_NSKeywords(i)
               success = .FALSE. 
            END IF  
         END DO  
         
      END SUBROUTINE CheckPhysicsInputIntegrity
!
!    **********       
     END MODULE PhysicsStorage_NS

