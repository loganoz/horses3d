!
!//////////////////////////////////////////////////////
!
!   @File:    PhysicsStorage_iNS.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Jun 19 17:39:26 2018
!   @Last revision date: Sat Jun 23 10:20:35 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: fce351220409e80ce5df1949249c2b870dd847aa
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
      Module Physics_iNSKeywordsModule
         IMPLICIT NONE 
         INTEGER, PARAMETER :: KEYWORD_LENGTH = 132
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: REYNOLDS_NUMBER_KEY       = "reynolds number"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: FROUDE_NUMBER_KEY         = "froude number"  
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: GRAVITY_DIRECTION_KEY     = "gravity direction"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: AOA_THETA_KEY             = "aoa theta"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: AOA_PHI_KEY               = "aoa phi"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: RIEMANN_SOLVER_NAME_KEY   = "riemann solver"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: LAMBDA_STABILIZATION_KEY  = "lambda stabilization"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: COMPUTE_GRADIENTS_KEY     = "compute gradients"
         
         CHARACTER(LEN=KEYWORD_LENGTH), DIMENSION(1) :: physics_iNSKeywords = [REYNOLDS_NUMBER_KEY]

         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: CENTRAL_SOLVER_NAME       = "central"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: LAXFRIEDRICHS_SOLVER_NAME = "lax-friedrichs"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: EXACT_SOLVER_NAME       = "exact"

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
      END MODULE Physics_iNSKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
!    
!    ******
     MODULE PhysicsStorage_iNS
!    ******
!
     USE SMConstants
     use FluidData_iNS
     use FileReadingUtilities, only: getArrayFromString
     
     IMPLICIT NONE

     private
     public    NINC
     public    INSRHO, INSU, INSV, INSW, INSP
     public    lambdaStab, computeGradients, whichRiemannSolver, whichAverage
     public    RIEMANN_CENTRAL, RIEMANN_LXF, RIEMANN_EXACT
     public    STANDARD_SPLIT
      
     public    ConstructPhysicsStorage_iNS, DestructPhysicsStorage_iNS, DescribePhysicsStorage_iNS
     public    CheckPhysics_iNSInputIntegrity
!
!    ----------------------------
!    Either NavierStokes or Euler
!    ----------------------------
!
     logical, protected :: computeGradients   = .true.
!
!    --------------------------
!!   The sizes of the NS system
!    --------------------------
!
     INTEGER, PARAMETER :: NINC = 5
!
!    -------------------------------------------
!!   The positions of the conservative variables
!    -------------------------------------------
!
     enum, bind(C) 
        enumerator :: INSRHO = 1, INSU, INSV, INSW, INSP
     end enum
!
!    --------------------------
!    Riemann solver definitions
!    --------------------------
!
     enum, bind(C)
        enumerator :: RIEMANN_CENTRAL = 1, RIEMANN_LXF, RIEMANN_EXACT
     end enum
     integer, protected :: whichRiemannSolver = -1
!
!    -----------------------------
!    Available averaging functions
!    -----------------------------
!
     enum, bind(C)
        enumerator :: STANDARD_SPLIT = 1
     end enum 
     integer            :: whichAverage               = -1
!
!    -------------------------------------
!    Lambda stabilization - 1.0 by default
!    -------------------------------------
!
     real(kind=RP), protected :: lambdaStab = 1.0_RP     

     type(Thermodynamics_t), target, private    :: ThermodynamicsAir = Thermodynamics_t("Air",1000.0_RP)
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
      SUBROUTINE ConstructPhysicsStorage_iNS( controlVariables, Lref, timeref, success )
      USE FTValueDictionaryClass
      USE Physics_iNSKeywordsModule
      use Utilities, only: toLower, almostEqual
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(FTValueDictionary)      :: controlVariables
      real(kind=RP), intent(inout) :: Lref, timeref
      LOGICAL                      :: success
!
!     ---------------
!     Local variables
!     ---------------
!
      CHARACTER(LEN=KEYWORD_LENGTH) :: keyword
      type(Thermodynamics_t), pointer  :: thermodynamics_
      type(RefValues_t)                :: refValues_
      type(Dimensionless_t)            :: dimensionless_
      real(kind=RP), allocatable       :: array(:)

!
!     --------------------
!     Collect input values
!     --------------------
!
      success = .TRUE.
      CALL CheckPhysics_iNSInputIntegrity(controlVariables,success)
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
!     *********************
!     Select flow equations
!     *********************
!
      dimensionless_ % Re = controlVariables % doublePrecisionValueForKey(REYNOLDS_NUMBER_KEY) 
!
!     ------------------------------------------------
!     Set molecular viscosity and thermal conductivity
!     ------------------------------------------------
!
      if ( .not. almostEqual(dimensionless_ % Re, 0.0_RP) ) then
         dimensionless_ % mu   = 1.0_RP / dimensionless_ % Re

      else
         dimensionless_ % mu = 0.0_RP

      end if
!
!     **************************************
!     Check if state gradients are requested
!     **************************************
!
      computeGradients = .true.
!
!     ********************
!     Set reference values: TODO read from parameter file
!                           Ok, but be sure to change the mesh reading accordingly (x = x / refValues % L)
!     ********************
!
      Lref = 1.0_RP
      refValues_ % rho = 1.0_RP
      refValues_ % V = 1.0_RP
      refValues_ % p = refValues_ % rho * POW2( refValues_ % V )
      refValues_ % mu = refValues_ % rho * refValues_ % V * Lref * dimensionless_ % mu
      timeref = Lref / refValues_ % V
!
!     ****************************************************
!     GRAVITY: Set the Froude number and gravity direction
!     ****************************************************
!
      if ( controlVariables % ContainsKey(FROUDE_NUMBER_KEY) ) then
         dimensionless_ % Fr = controlVariables % DoublePrecisionValueForKey(FROUDE_NUMBER_KEY)

      else
!
!        Default Froude number: earth's gravity
!        --------------------------------------
         dimensionless_ % Fr = refValues_ % V / sqrt(9.81_RP * Lref)

      end if

      if ( controlVariables % ContainsKey(GRAVITY_DIRECTION_KEY) ) then
         array = GetArrayFromString( controlVariables % StringValueForKey(GRAVITY_DIRECTION_KEY,&
                                                                             KEYWORD_LENGTH))
         dimensionless_ % gravity_dir = array(1:3)

         if ( norm2(dimensionless_ % gravity_dir) < epsilon(1.0_RP)*10.0_RP ) then
!
!           Disable gravity
!           ---------------
            dimensionless_ % gravity_dir = 0.0_RP
            dimensionless_ % invFroudeSquare = 0.0_RP
            dimensionless_ % Fr = huge(1.0_RP)
         else
            dimensionless_ % gravity_dir = dimensionless_ % gravity_dir / norm2(dimensionless_ % gravity_dir)
            dimensionless_ % invFroudeSquare = 1.0_RP / POW2(dimensionless_ % Fr)

         end if
      else
         if ( controlVariables % ContainsKey(FROUDE_NUMBER_KEY) ) then
            print*, "When specifying a Froude number, the gravity direction must be specified"
            print*, "Gravity direction = [x,y,z]"
            errorMessage(STD_OUT)
            stop

         else
!
!           Gravity is disabled
!           -------------------
            dimensionless_ % gravity_dir = 0.0_RP
            dimensionless_ % Fr = huge(1.0_RP)
            dimensionless_ % invFroudeSquare = 0.0_RP

         end if
      end if
!
!     *********************************************
!     Choose the Riemann solver (by default is ERS)
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
         case(CENTRAL_SOLVER_NAME) 
            whichRiemannSolver = RIEMANN_CENTRAL

         case(LAXFRIEDRICHS_SOLVER_NAME)
            whichRiemannSolver = RIEMANN_LXF 

         case(EXACT_SOLVER_NAME)
            whichRiemannSolver = RIEMANN_EXACT
   
         case default 
            print*, "Riemann solver: ", trim(keyword), " is not implemented."
            print*, "Options available are:"
            print*, "   * Central"
            print*, "   * Lax-Friedrichs"
            print*, "   * Exact"
            errorMessage(STD_OUT)
            stop
         end select 
      else 
!
!        Select Roe by default
!        ---------------------
         whichRiemannSolver = RIEMANN_EXACT 

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
      CALL DescribePhysicsStorage_iNS()

      END SUBROUTINE ConstructPhysicsStorage_iNS
!
!     ///////////////////////////////////////////////////////
!
!     -------------------------------------------------
!!    Destructor: Does nothing for this storage
!     -------------------------------------------------
!
      SUBROUTINE DestructPhysicsStorage_iNS
      
      END SUBROUTINE DestructPhysicsStorage_iNS
!
!     //////////////////////////////////////////////////////
!
!     -----------------------------------------
!!    Descriptor: Shows the gathered data
!     -----------------------------------------
!
      SUBROUTINE DescribePhysicsStorage_iNS()
         USE Headers
         use MPI_Process_Info
         IMPLICIT NONE
         real(kind=RP)  :: pRef

         if ( .not. MPI_Process % isRoot ) return 

         write(STD_OUT,'(/,/)')
         call Section_Header("Loading incompressible Navier-Stokes physics")

         write(STD_OUT,'(/)')
         call SubSection_Header("Fluid data")
         write(STD_OUT,'(30X,A,A22,A10)') "->" , "Fluid: " , "Air"

         write(STD_OUT,'(/)')
         call SubSection_Header("Reference quantities")
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference pressure: " , refValues % p, " Pa."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference density: " , refValues % rho , " kg/m^3."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference velocity: " , refValues % V , " m/s."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference viscosity: ",refValues % mu , " Pa·s."

         write(STD_OUT,'(/)')
         call SubSection_Header("Dimensionless quantities")
         write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Reynolds number: " , dimensionless % Re
         write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Froude number: " , dimensionless % Fr
         write(STD_OUT,'(30X,A,A20,A,F4.1,A,F4.1,A,F4.1,A)') "->" , "Gravity direction: ","[", &
                                                   dimensionless % gravity_dir(1), ", ", &
                                                   dimensionless % gravity_dir(2), ", ", &
                                                   dimensionless % gravity_dir(3), "]"

      END SUBROUTINE DescribePhysicsStorage_iNS
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE CheckPhysics_iNSInputIntegrity( controlVariables, success )  
         USE FTValueDictionaryClass
         USE Physics_iNSKeywordsModule
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
         
         DO i = 1, SIZE(physics_iNSKeywords)
            obj => controlVariables % objectForKey(physics_iNSKeywords(i))
            IF ( .NOT. ASSOCIATED(obj) )     THEN
               PRINT *, "Input file is missing entry for keyword: ",physics_iNSKeywords(i)
               success = .FALSE. 
            END IF  
         END DO  
         
      END SUBROUTINE CheckPhysics_iNSInputIntegrity
!
!    **********       
     END MODULE PhysicsStorage_iNS
!    **********

