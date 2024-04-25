#include "Includes.h"
      Module Physics_iNSKeywordsModule
         IMPLICIT NONE
         INTEGER, PARAMETER :: KEYWORD_LENGTH = 132
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: REFERENCE_VELOCITY_KEY         = "reference velocity (m/s)"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: NUMBER_OF_FLUIDS_KEY           = "number of fluids (1/2)"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: MAXIMUM_DENSITY_KEY            = "maximum density (kg/m^3)"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: MINIMUM_DENSITY_KEY            = "minimum density (kg/m^3)"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: ARTIFICIAL_COMPRESSIBILITY_KEY = "artificial compressibility factor"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: GRAVITY_ACCELERATION_KEY       = "gravity acceleration (m/s^2)"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: GRAVITY_DIRECTION_KEY          = "gravity direction"
!
!        *****************
!        Mode with 1 fluid
!        *****************
!
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: DENSITY_KEY    =  "density (kg/m^3)"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: VISCOSITY_KEY  =  "viscosity (pa.s)"
!
!        *****************
!        Mode with 2 fluid
!        *****************
!
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: FLUID1_DENSITY_KEY    =  "fluid 1 density (kg/m^3)"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: FLUID2_DENSITY_KEY    =  "fluid 2 density (kg/m^3)"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: FLUID1_VISCOSITY_KEY  =  "fluid 1 viscosity (pa.s)"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: FLUID2_VISCOSITY_KEY  =  "fluid 2 viscosity (pa.s)"

         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: AOA_THETA_KEY             = "aoa theta"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: AOA_PHI_KEY               = "aoa phi"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: COMPUTE_GRADIENTS_KEY     = "compute gradients"

         ! !PARTICLES 
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: particlesKey             = "lagrangian particles"         
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: numberOfParticlesKey     = "number of particles" 
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: particlesPerParcelKey    = "particles per parcel"     
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: sourceTermKey            = "high order particles source term"     
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: STOKES_NUMBER_PART_KEY   = "stokes number" 
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: GAMMA_PART_KEY           = "gamma" 
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PHI_M_PART_KEY           = "phi_m" 
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: I0_PART_KEY              = "radiation source" 
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: MIN_BOX_KEY              = "minimum box" 
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: MAX_BOX_KEY              = "maximum box" 
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: BC_BOX_KEY               = "bc box" 
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PART_FILE_KEY            = "particles file"
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PART_LOG_FILE_KEY        = "vel and temp from file"
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PART_LOG_INJ_KEY         = "injection"         
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PART_INJ_KEY             = "particles injection"
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PART_NUMB_PER_STEP_KEY   = "particles per step"
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: PART_PERIOD_KEY          = "particles iter period"
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: INJ_VEL_KEY              = "particles injection velocity"
         ! CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: INJ_TEMP_KEY             = "particles injection temperature"     
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
     use FileReadingUtilities, only: getRealArrayFromString

     IMPLICIT NONE

     private
     public    NCONS, NGRAD
     public    INSRHO, INSRHOU, INSRHOV, INSRHOW, INSP
     public    computeGradients
     public    enableGravity

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
     INTEGER, PARAMETER :: NCONS = 5, NGRAD = 5
!
!    -------------------------------------------
!!   The positions of the conservative variables
!    -------------------------------------------
!
     enum, bind(C)
        enumerator :: INSRHO = 1, INSRHOU, INSRHOV, INSRHOW, INSP
     end enum

     logical, protected       :: enableGravity         = .false.
!
!    ========
     contains
!    ========
!
!     ///////////////////////////////////////////////////////
!
!     --------------------------------------------------
!     Constructor: Define default values for the physics
!     variables.
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
      real(kind=RP), intent(in)    :: Lref
      real(kind=RP), intent(out)   :: timeref
      LOGICAL                      :: success
!
!     ---------------
!     Local variables
!     ---------------
!
      CHARACTER(LEN=KEYWORD_LENGTH)   :: keyword
      type(Thermodynamics_t)          :: thermodynamics_

      type(RefValues_t)               :: refValues_
      type(Dimensionless_t)           :: dimensionless_
      real(kind=RP)                   :: array(3)
!
!     --------------------
!     Collect input values
!     --------------------
!
      success = .TRUE.
      CALL CheckPhysics_iNSInputIntegrity(controlVariables,success)
      IF(.NOT. success) RETURN
!
!     **************************************
!     Check if state gradients are requested
!     **************************************
!
      computeGradients = .true.
!
!     ******************
!     Set thermodynamics
!     ******************
!
      thermodynamics_ % number_of_fluids = controlVariables % IntegerValueForKey(NUMBER_OF_FLUIDS_KEY)

      allocate(thermodynamics_ % rho(thermodynamics_ % number_of_fluids))
      allocate(thermodynamics_ % mu (thermodynamics_ % number_of_fluids))

      select case(thermodynamics_ % number_of_fluids)
      case(1)
         thermodynamics_ % rho(1) = controlVariables % DoublePrecisionValueForKey(DENSITY_KEY)
         thermodynamics_ % mu (1) = controlVariables % DoublePrecisionValueForKey(VISCOSITY_KEY)

      case(2)
         thermodynamics_ % rho(1) = controlVariables % DoublePrecisionValueForKey(FLUID1_DENSITY_KEY)
         thermodynamics_ % rho(2) = controlVariables % DoublePrecisionValueForKey(FLUID2_DENSITY_KEY)

         thermodynamics_ % mu (1) = controlVariables % DoublePrecisionValueForKey(FLUID1_VISCOSITY_KEY)
         thermodynamics_ % mu (2) = controlVariables % DoublePrecisionValueForKey(FLUID2_VISCOSITY_KEY)

      end select


      if ( controlVariables % ContainsKey(MAXIMUM_DENSITY_KEY) ) then
         thermodynamics_ % rho_max = controlVariables % DoublePrecisionValueForKey(MAXIMUM_DENSITY_KEY)
      else
         thermodynamics_ % rho_max = huge(1.0_RP)

      end if

      if ( controlVariables % ContainsKey(MINIMUM_DENSITY_KEY) ) then
         thermodynamics_ % rho_min = controlVariables % DoublePrecisionValueForKey(MINIMUM_DENSITY_KEY)
      else
         thermodynamics_ % rho_min = -huge(1.0_RP)

      end if

!
!     ********************
!     Set reference values
!     ********************
!
      refValues_ % rho = thermodynamics_ % rho(1)
      refValues_ % V   = controlVariables % DoublePrecisionValueForKey(REFERENCE_VELOCITY_KEY)
      refValues_ % p   = refValues_ % rho * POW2( refValues_ % V )
      refValues_ % mu  = thermodynamics_ % mu(1)
      timeref          = Lref / refValues_ % V
!
!     ****************************
!     Set dimensionless quantities
!     ****************************
!
      allocate(dimensionless_ % rho(thermodynamics_ % number_of_fluids))
      allocate(dimensionless_ % mu (thermodynamics_ % number_of_fluids))

      dimensionless_ % rho = thermodynamics_ % rho / refValues_ % rho
      dimensionless_ % mu(1) = thermodynamics_ % mu(1) / (thermodynamics_ % rho(1) * refValues_ % V * Lref)

      if ( thermodynamics_ % number_of_fluids .eq. 2 ) then
         if ( .not. almostEqual(thermodynamics_ % mu(1), 0.0_RP)) then
            dimensionless_ % mu(2) = dimensionless_ % mu(1) * thermodynamics_ % mu(2) / thermodynamics_ % mu(1)
            dimensionless_ % Re = 1.0_RP / max ( dimensionless_ % mu(1), epsilon(1.0_RP) )
         else
            dimensionless_ % mu(2) = 0.0_RP
            dimensionless_ % Re    = 0.0_RP

         end if
      else
         dimensionless_ % Re = 1.0_RP / max ( dimensionless_ % mu(1), epsilon(1.0_RP) )
      end if
!
!     **************************
!     Artificial compressibility
!     **************************
!
      thermodynamics_ % rho0c02 = maxval(dimensionless_ % rho) * controlVariables % DoublePrecisionValueForKey(ARTIFICIAL_COMPRESSIBILITY_KEY)
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
!     *************
!     Gravity force
!     *************
!
      if ( controlVariables % ContainsKey(GRAVITY_DIRECTION_KEY) ) then
         array = getRealArrayFromString( controlVariables % StringValueForKey(GRAVITY_DIRECTION_KEY,&
                                                                             KEYWORD_LENGTH))
         dimensionless_ % gravity_dir = array(1:3) / norm2(array(1:3))
      end if

      refValues_ % g0 = controlVariables % DoublePrecisionValueForKey(GRAVITY_ACCELERATION_KEY)

      if ( almostEqual(abs(refValues_ % g0), 0.0_RP) ) then
         enableGravity = .false.
         dimensionless_ % Fr = 0.0_RP
         dimensionless_ % invFr2 = 0.0_RP

      else
         enableGravity = .true.
         dimensionless_ % invFr2 = refValues_ % g0 * Lref / POW2(refValues_ % V)
         dimensionless_ % Fr     = 1.0_RP / sqrt(dimensionless_ % invFr2)

      end if
!
!     **********************************************************************
!     Set the global (proteted) thermodynamics, dimensionless, and refValues
!     **********************************************************************
!
      call setThermodynamics(thermodynamics_)
      call setDimensionless (dimensionless_ )
      call setRefValues     (refValues_     )

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
         write(STD_OUT,'(30X,A,A22,I0)') "->" , "Number of fluids: " , thermodynamics % number_of_fluids
         select case(thermodynamics % number_of_fluids)
         case(1)
            write(STD_OUT,'(30X,A,A22,F10.3,A)') "->" , "Fluid density: " , thermodynamics % rho(1), " kg/m^3"
            write(STD_OUT,'(30X,A,A22,F10.3,A)') "->" , "Fluid viscosity: " , thermodynamics % mu(1), " Pa.s"

         case(2)
            write(STD_OUT,'(30X,A,A22,F10.3,A)') "->" , "Fluid 1 density: " , thermodynamics % rho(1), " kg/m^3"
            write(STD_OUT,'(30X,A,A22,F10.3,A)') "->" , "Fluid 2 density: " , thermodynamics % rho(2), " kg/m^3"
            write(STD_OUT,'(30X,A,A22,1pG10.3,A)') "->" , "Fluid 1 viscosity: " , thermodynamics % mu(1), " Pa.s"
            write(STD_OUT,'(30X,A,A22,1pG10.3,A)') "->" , "Fluid 2 viscosity: " , thermodynamics % mu(2), " Pa.s"
         end select

         write(STD_OUT,'(30X,A,A22,F10.3,A)') "->" , "Artificial compressibility " , thermodynamics % rho0c02, "-"

         write(STD_OUT,'(/)')
         call SubSection_Header("Reference quantities")
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference pressure: " , refValues % p, " Pa."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference density: " , refValues % rho , " kg/m^3."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference velocity: " , refValues % V , " m/s."
         write(STD_OUT,'(30X,A,A30,1pG10.3,A)') "->" , "Reference viscosity: ",refValues % mu , " Pa·s."

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
!
!        *******************************************************************
!           In this solver there are not compulsory keywords, but they
!           are given default values
!        *******************************************************************
!
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
         INTEGER                  :: i, nF
         real(kind=RP)            :: array(3)
         success = .TRUE.
         
!         DO i = 1, SIZE(physics_iNSKeywords)
!            obj => controlVariables % objectForKey(physics_iNSKeywords(i))
!            IF ( .NOT. ASSOCIATED(obj) )     THEN
!               PRINT *, "Input file is missing entry for keyword: ",physics_iNSKeywords(i)
!               success = .FALSE. 
!            END IF  
!         END DO  

         if ( .not. controlVariables % ContainsKey(REFERENCE_VELOCITY_KEY) ) then
            call controlVariables % AddValueForKey("1.0", REFERENCE_VELOCITY_KEY)
         end if

         if ( .not. controlVariables % ContainsKey(NUMBER_OF_FLUIDS_KEY) ) then
            call controlVariables % AddValueForKey("1", NUMBER_OF_FLUIDS_KEY)
         end if

         if ( .not. controlVariables % ContainsKey(ARTIFICIAL_COMPRESSIBILITY_KEY) ) then
            call controlVariables % AddValueForKey("1000.0", ARTIFICIAL_COMPRESSIBILITY_KEY)
         end if

         nF = controlVariables % IntegerValueForKey(NUMBER_OF_FLUIDS_KEY)

         select case (nF)
         case (1)
            if ( .not. controlVariables % ContainsKey(DENSITY_KEY)) then
               call controlVariables % AddValueForKey("1.0", DENSITY_KEY)
            end if

            if ( .not. controlVariables % ContainsKey(VISCOSITY_KEY)) then
               call controlVariables % AddValueForKey("0.0", VISCOSITY_KEY)
            end if
         case (2)
            if ( .not. controlVariables % ContainsKey(FLUID1_DENSITY_KEY)) then
               print*, "Specify density for fluid #1 using:"
               print*, "   ",trim(FLUID1_DENSITY_KEY), " = #value"
               errorMessage(STD_OUT)
               error stop
            end if

            if ( .not. controlVariables % ContainsKey(FLUID2_DENSITY_KEY)) then
               print*, "Specify density for fluid #2 using:"
               print*, "   ",trim(FLUID2_DENSITY_KEY), " = #value"
               errorMessage(STD_OUT)
               error stop
            end if

            if ( .not. controlVariables % ContainsKey(FLUID1_VISCOSITY_KEY)) then
               call controlVariables % AddValueForKey("0.0", FLUID1_VISCOSITY_KEY)
            end if

            if ( .not. controlVariables % ContainsKey(FLUID2_VISCOSITY_KEY)) then
               call controlVariables % AddValueForKey("0.0", FLUID2_VISCOSITY_KEY)
            end if

         end select
!
!        *************
!        Gravity force
!        *************
!
         if ( controlVariables % ContainsKey(GRAVITY_DIRECTION_KEY) ) then
            array = getRealArrayFromString( controlVariables % StringValueForKey(GRAVITY_DIRECTION_KEY,&
                                                                                KEYWORD_LENGTH))
            if ( norm2(array) < epsilon(1.0_RP)*10.0_RP ) then
!
!              Error
!              -----
               print*, "Incorrect gravity direction vector"
               errorMessage(STD_OUT)
               error stop

            end if

            if ( .not. controlVariables % ContainsKey(GRAVITY_ACCELERATION_KEY) ) then
               call controlVariables % AddValueForKey("9.81d0", GRAVITY_ACCELERATION_KEY)

            end if

         else
            if ( controlVariables % ContainsKey(GRAVITY_ACCELERATION_KEY) ) then
               print*, "Gravity acceleration requires gravity direction."
               print*, "Specify gravity direction with:"
               print*, "     ", GRAVITY_DIRECTION_KEY, " = [x,y,z]"
               errorMessage(STD_OUT)
               error stop

            else
               call controlVariables % AddValueForKey("0.0d0", GRAVITY_ACCELERATION_KEY)

            end if

         end if


      END SUBROUTINE CheckPhysics_iNSInputIntegrity
!
!    **********
     END MODULE PhysicsStorage_iNS
!    **********
