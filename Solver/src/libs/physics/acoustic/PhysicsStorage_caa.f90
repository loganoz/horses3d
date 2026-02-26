#include "Includes.h"
      Module Physics_CAAKeywordsModule
         IMPLICIT NONE
         INTEGER, PARAMETER :: KEYWORD_LENGTH = 132
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: REFERENCE_TEMPERATURE_KEY      = "reference temperature (k)"
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: REFERENCE_PRESSURE_KEY         = "reference pressure (pa)"
         character(LEN = KEYWORD_LENGTH), parameter :: REFERENCE_VELOCITY_KEY         = "reference velocity (m/s)"
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: MACH_NUMBER_KEY                = "mach number"
         CHARACTER(LEN = KEYWORD_LENGTH), PARAMETER :: FLOW_EQUATIONS_KEY             = "flow equations"

         character(len=KEYWORD_LENGTH), parameter :: ARTIFICIAL_COMPRESSIBILITY_KEY = "artificial sound speed square (m/s)"
         character(len=KEYWORD_LENGTH), parameter :: FLUID1_COMPRESSIBILITY_KEY = "fluid 1 sound speed square (m/s)"
         character(len=KEYWORD_LENGTH), parameter :: FLUID2_COMPRESSIBILITY_KEY = "fluid 2 sound speed square (m/s)"

         CHARACTER(LEN=KEYWORD_LENGTH), DIMENSION(2) :: physics_CAAKeywords = [MACH_NUMBER_KEY, FLOW_EQUATIONS_KEY]

      END MODULE Physics_CAAKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
!
!    ******
     MODULE PhysicsStorage_CAA
!    ******
!
     USE SMConstants
     use FluidData_CAA
     use FileReadingUtilities, only: getRealArrayFromString

     IMPLICIT NONE

     private
     public    NCONS, NGRAD, NCONSB
     public    ICAARHO, ICAAU, ICAAV, ICAAW, ICAAP
     public    IBRHO, IBU, IBV, IBW, IBP, IBA2
     public    IRHO, IRHOU, IRHOV, IRHOW, IRHOE ! NS solver
     public    INSRHO, INSRHOU, INSRHOV, INSRHOW, INSP ! iNS solver
     public    IMC, IMSQRHOU, IMSQRHOV, IMSQRHOW, IMP ! MU solver
     public    computeGradients, flowIsNavierStokes

     public    ConstructPhysicsStorage_CAA, DestructPhysicsStorage_CAA, DescribePhysicsStorage_CAA
     public    CheckPhysicsCAAInputIntegrity
     public    GRADVARS_STATE, grad_vars
!
!    ----------------------------
!    Either NavierStokes or Euler, for compatibility reasons
!    ----------------------------
!
     logical, protected :: flowIsNavierStokes = .false.
     logical, protected :: computeGradients   = .false.
!
!    --------------------------
!!   The sizes of the CAA system
!    --------------------------
!
     INTEGER, PARAMETER :: NCONS = 5, NGRAD = 5
     ! for APE and LEE we have one more variable, the speed of sound squared
     INTEGER, PARAMETER ::  NCONSB = NCONS + 1
!
!    -------------------------------------------
!!   The positions of the conserv variables
!    -------------------------------------------
!
     INTEGER, PARAMETER       :: ICAARHO = 1 , ICAAU = 2 , ICAAV = 3 , ICAAW = 4 , ICAAP = 5
!
!!   The positions of the primitive variables of the base flow
     INTEGER, PARAMETER       :: IBRHO = 1 , IBU = 2 , IBV = 3 , IBW = 4 , IBP = 5, IBA2 = 6

!!   The positions of the conservative variables of the base flow from the NS solver
     INTEGER, PARAMETER       :: IRHO = 1 , IRHOU = 2 , IRHOV = 3 , IRHOW = 4 , IRHOE = 5

!!   The positions of the conservative variables of the base flow from the iNS solver
     enum, bind(C)
        enumerator :: INSRHO = 1, INSRHOU, INSRHOV, INSRHOW, INSP
     end enum

!!   The positions of the conservative variables of the base flow from the MU solver
     INTEGER, PARAMETER       :: IMC = 1 , IMSQRHOU = 2 , IMSQRHOV = 3 , IMSQRHOW = 4 , IMP = 5
!
!    --------------------------------
!    Choice of the gradient variables
!    --------------------------------
!
     enum, bind(C)
        enumerator :: GRADVARS_STATE
     end enum
     integer, protected :: grad_vars = GRADVARS_STATE

!
!    ------------------------------
!    Thermodynamic specs of the Air
!    ------------------------------
!
     type(Thermodynamics_t), target, private :: ThermodynamicsAir = Thermodynamics_t( &
                                                              "Air", & ! Name
                                                      287.15_RP, & ! R  ! J / (kg * °K)
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
                       287.15_RP * 1.4_RP / ( 1.4_RP - 1.0_RP ), & ! cp
                               287.15_RP / ( 1.4_RP - 1.0_RP ),  & ! cv
                                 [ 117649.0_RP, 117649.0_RP ]    & ! c02
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
      SUBROUTINE ConstructPhysicsStorage_CAA( controlVariables, Lref, timeref, success )
      USE FTValueDictionaryClass
      USE Physics_CAAKeywordsModule
      use Utilities, only: toLower, almostEqual
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(FTValueDictionary)      :: controlVariables
      real(kind=RP), intent(in)    :: Lref
      real(kind=RP), intent(out)   :: timeref
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
      CALL CheckPhysicsCAAInputIntegrity(controlVariables,success)
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
!
!     *********************
!     Select flow equations
!     *********************
!
      keyword = controlVariables % stringValueForKey(FLOW_EQUATIONS_KEY,KEYWORD_LENGTH)
      CALL toLower(keyword)

      ! if ( keyword == "ape" ) then
      ! end if
!
!     **************************************************
!     Set reference values with the base flow variables
!     **************************************************
!
      if ( controlVariables % ContainsKey(REFERENCE_TEMPERATURE_KEY) ) then
         refValues_ % T = controlVariables % DoublePrecisionValueForKey(REFERENCE_TEMPERATURE_KEY)
      else
         refValues_ % T = 520.0_RP*5.0_RP/9.0_RP     ! ≃ 288.88 K
      end if

      if ( controlVariables % ContainsKey(REFERENCE_PRESSURE_KEY) ) then
         refValues_ % rho = controlVariables % DoublePrecisionValueForKey(REFERENCE_PRESSURE_KEY) / (thermodynamics_ % R * refValues_ % T)
      else
         refValues_ % rho = 101325.0_RP / (thermodynamics_ % R * refValues_ % T)
      end if

      if ( controlVariables % ContainsKey(MACH_NUMBER_KEY) ) then
            dimensionless_ % Mach = controlVariables % doublePrecisionValueForKey(MACH_NUMBER_KEY)
            refValues_ % V =   dimensionless_ % Mach &
                        * sqrt( thermodynamics_ % gamma * thermodynamics_ % R * refValues_ % T )
      elseif ( controlVariables % ContainsKey(REFERENCE_VELOCITY_KEY) ) then
            refValues_ % V = controlVariables % doublePrecisionValueForKey(REFERENCE_VELOCITY_KEY)
            dimensionless_ % Mach = refValues_ % V &
                        / sqrt( thermodynamics_ % gamma * thermodynamics_ % R * refValues_ % T )
      else
            print*, "*** ERROR: Specify Mach number or reference velocity"
            error stop
      end if

      
      if ( controlVariables % ContainsKey(ARTIFICIAL_COMPRESSIBILITY_KEY) ) then
         thermodynamics_ % c02 = controlVariables % DoublePrecisionValueForKey(ARTIFICIAL_COMPRESSIBILITY_KEY)
      else
!      - Default to 1000.0
         thermodynamics_ % c02 = 1000.0_RP

         if ( controlVariables % ContainsKey(FLUID1_COMPRESSIBILITY_KEY) ) then
            thermodynamics_ % c02(1) = controlVariables % DoublePrecisionValueForKey(FLUID1_COMPRESSIBILITY_KEY)
         end if

         if ( controlVariables % ContainsKey(FLUID2_COMPRESSIBILITY_KEY) ) then
            thermodynamics_ % c02(2) = controlVariables % DoublePrecisionValueForKey(FLUID2_COMPRESSIBILITY_KEY)
         end if
         
      end if

      refValues_ % p = refValues_ % rho * POW2( refValues_ % V )

      timeref = Lref / refValues_ % V

      dimensionless_ % gammaM2 = thermodynamics_ % gamma * POW2( dimensionless_ % Mach )
!
!     **********************************************************************
!     Set the global (proteted) thermodynamics, dimensionless, and refValues
!     **********************************************************************
!
      call setThermodynamics( thermodynamics_ )
      call setDimensionless( dimensionless_ )
      call setRefValues( refValues_ )

      END SUBROUTINE ConstructPhysicsStorage_CAA
!
!     ///////////////////////////////////////////////////////
!
!     -------------------------------------------------
!!    Destructor: Does nothing for this storage
!     -------------------------------------------------
!
      SUBROUTINE DestructPhysicsStorage_CAA

      END SUBROUTINE DestructPhysicsStorage_CAA
!
!     //////////////////////////////////////////////////////
!
!     -----------------------------------------
!!    Descriptor: Shows the gathered data
!     -----------------------------------------
!
      SUBROUTINE DescribePhysicsStorage_CAA()
         USE Headers
         use MPI_Process_Info
         IMPLICIT NONE
         real(kind=RP)  :: pRef

         if ( .not. MPI_Process % isRoot ) return

         pRef = thermodynamics % R * refValues % rho * refValues % T

         write(STD_OUT,'(/,/)')
         call Section_Header("Loading APE physics")

         write(STD_OUT,'(/)')
         call SubSection_Header("Fluid data")
         write(STD_OUT,'(30X,A,A22,A10)') "->" , "Gas: " , "Air"
         write(STD_OUT,'(30X,A,A22,F10.3,A)') "->" , "State constant: " , thermodynamics % R, " I.S."
         write(STD_OUT,'(30X,A,A22,F10.3)') "->" , "Specific heat ratio: " , thermodynamics % gamma

         write(STD_OUT,'(/)')
         call SubSection_Header("Reference quantities of base flow")
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference Temperature: " , refValues % T, " K."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference pressure: " , pRef, " Pa."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference density: " , refValues % rho , " kg/m^3."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference velocity: " , refValues % V , " m/s."

         write(STD_OUT,'(/)')
         call SubSection_Header("Dimensionless quantities of base flow")
         write(STD_OUT,'(30X,A,A27,F10.3)') "->" , "Mach number: " , dimensionless % Mach

      END SUBROUTINE DescribePhysicsStorage_CAA
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE CheckPhysicsCAAInputIntegrity( controlVariables, success )
         USE FTValueDictionaryClass
         USE Physics_CAAKeywordsModule
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

         DO i = 1, SIZE(physics_CAAKeywords)
            obj => controlVariables % objectForKey(physics_CAAKeywords(i))
            IF ( .NOT. ASSOCIATED(obj) )     THEN
               PRINT *, "Input file is missing entry for keyword: ",physics_CAAKeywords(i)
               success = .FALSE.
            END IF
         END DO

      END SUBROUTINE CheckPhysicsCAAInputIntegrity

      subroutine SetGradientVariables(grad_vars_)
         implicit none
         integer, intent(in)  :: grad_vars_

         select case(grad_vars_)
         case(GRADVARS_STATE )
            grad_vars = grad_vars_
         case default
            print*, "Unrecognized option"
            errorMessage(STD_OUT)
            error stop
         end select

      end subroutine SetGradientVariables

!
!    **********
     END MODULE PhysicsStorage_CAA
!    **********
