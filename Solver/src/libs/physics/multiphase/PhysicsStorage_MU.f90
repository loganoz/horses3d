#include "Includes.h"
      Module Physics_MUKeywordsModule
         IMPLICIT NONE
         INTEGER, parameter :: KEYWORD_LENGTH = 132
         character(len=KEYWORD_LENGTH), parameter :: REFERENCE_VELOCITY_KEY         = "reference velocity (m/s)"
         character(len=KEYWORD_LENGTH), parameter :: MAXIMUM_DENSITY_KEY            = "maximum density (kg/m^3)"
         character(len=KEYWORD_LENGTH), parameter :: MINIMUM_DENSITY_KEY            = "minimum density (kg/m^3)"
         character(len=KEYWORD_LENGTH), parameter :: ARTIFICIAL_COMPRESSIBILITY_KEY = "artificial sound speed square (m/s)"
         character(len=KEYWORD_LENGTH), parameter :: GRAVITY_ACCELERATION_KEY       = "gravity acceleration (m/s^2)"
         character(len=KEYWORD_LENGTH), parameter :: GRAVITY_DIRECTION_KEY          = "gravity direction"
         character(len=KEYWORD_LENGTH), parameter :: VELOCITY_DIRECTION_KEY         = "velocity direction"

         character(len=KEYWORD_LENGTH), parameter :: FLUID1_DENSITY_KEY    =  "fluid 1 density (kg/m^3)"
         character(len=KEYWORD_LENGTH), parameter :: FLUID2_DENSITY_KEY    =  "fluid 2 density (kg/m^3)"
         character(len=KEYWORD_LENGTH), parameter :: FLUID1_VISCOSITY_KEY  =  "fluid 1 viscosity (pa.s)"
         character(len=KEYWORD_LENGTH), parameter :: FLUID2_VISCOSITY_KEY  =  "fluid 2 viscosity (pa.s)"

      END MODULE Physics_MUKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
!
!    ******
     MODULE PhysicsStorage_MU
!    ******
!
     USE SMConstants
     use FluidData_MU
     use FileReadingUtilities, only: getRealArrayFromString
     use Utilities,            only: toLower, almostEqual

     IMPLICIT NONE

     private
     public    NCONS, NGRAD
     public    IMC, IMSQRHOU, IMSQRHOV, IMSQRHOW, IMP
     public    IGMU, IGU, IGV, IGW, IGP
     public    computeGradients
     public    enableGravity

     public    ConstructPhysicsStorage_MU, DestructPhysicsStorage_MU, DescribePhysicsStorage_MU
     public    CheckPhysics_MUInputIntegrity
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
     INTEGER, parameter :: NCONS = 5, NGRAD = 5
!
!    -------------------------------------------
!!   The positions of the conservative variables
!    -------------------------------------------
!
     enum, bind(C)
        enumerator :: IMC = 1, IMSQRHOU, IMSQRHOV, IMSQRHOW, IMP
     end enum

     enum, bind(C)
        enumerator :: IGMU = 1, IGU, IGV, IGW, IGP
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
      SUBROUTINE ConstructPhysicsStorage_MU( controlVariables, Lref, timeref, success )
      USE FTValueDictionaryClass
      USE Physics_MUKeywordsModule
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
      character(len=KEYWORD_LENGTH)   :: keyword
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
      CALL CheckPhysics_MUInputIntegrity(controlVariables,success)
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
      if ( controlVariables % ContainsKey(FLUID1_DENSITY_KEY) ) then
         thermodynamics_ % rho(1) = controlVariables % DoublePrecisionValueForKey(FLUID1_DENSITY_KEY)
      else
!      - Default to 1.0
         thermodynamics_ % rho(1) = 1.0_RP
      end if


      if ( controlVariables % ContainsKey(FLUID2_DENSITY_KEY) ) then
         thermodynamics_ % rho(2) = controlVariables % DoublePrecisionValueForKey(FLUID2_DENSITY_KEY)
      else
!      - Default to 1.0
         thermodynamics_ % rho(2) = 1.0_RP
      end if

      if ( controlVariables % ContainsKey(FLUID1_VISCOSITY_KEY) ) then
         thermodynamics_ % mu(1) = controlVariables % DoublePrecisionValueForKey(FLUID1_VISCOSITY_KEY)
      else
!      - Default to 0.0
         thermodynamics_ % mu(1) = 0.0_RP
      end if

      if ( controlVariables % ContainsKey(FLUID2_VISCOSITY_KEY) ) then
         thermodynamics_ % mu(2) = controlVariables % DoublePrecisionValueForKey(FLUID2_VISCOSITY_KEY)
      else
!      - Default to 0.0
         thermodynamics_ % mu(2) = 0.0_RP
      end if

      if ( controlVariables % ContainsKey(ARTIFICIAL_COMPRESSIBILITY_KEY) ) then
         thermodynamics_ % c02 = controlVariables % DoublePrecisionValueForKey(ARTIFICIAL_COMPRESSIBILITY_KEY)
      else
!      - Default to 1000.0
         thermodynamics_ % c02 = 1000.0_RP
      end if
!
!     ********************
!     Set reference values
!     ********************
!
      refValues_ % rho = thermodynamics_ % rho(1)
      refValues_ % V   = controlVariables % DoublePrecisionValueForKey(REFERENCE_VELOCITY_KEY)
      refValues_ % p   = refValues_ % rho * POW2( refValues_ % V )

      if ( .not. almostEqual(thermodynamics_ % mu(1), 0.0_RP) ) then
         refValues_ % mu = thermodynamics_ % mu(1)
      else
         refValues_ % mu = thermodynamics_ % mu(2)
      end if

      refValues_ % g0  = controlVariables % DoublePrecisionValueForKey(GRAVITY_ACCELERATION_KEY)
      timeref          = Lref / refValues_ % V
!
!     ****************************
!     Set dimensionless quantities
!     ****************************
!
!     -------------
!     Thermodynamic
!     -------------
!
      dimensionless_ % rho = thermodynamics_ % rho / refValues_ % rho
      dimensionless_ % mu  = thermodynamics_ % mu  / (refValues_ % rho * refValues_ % V * Lref)

      if ( .not. almostEqual(dimensionless_ % mu(1), 0.0_RP) ) then
         dimensionless_ % Re(1) = 1.0_RP / dimensionless_ % mu(1)
      else
         dimensionless_ % Re(1) = 0.0_RP
      end if

      if ( .not. almostEqual(dimensionless_ % mu(2), 0.0_RP) ) then
         dimensionless_ % Re(2) = (dimensionless_ % rho(2)/dimensionless_ % rho(1)) / dimensionless_ % mu(2)
      else
         dimensionless_ % Re(2) = 0.0_RP
      end if
!
!     -------
!     Gravity
!     -------
!
      if ( controlVariables % ContainsKey(GRAVITY_DIRECTION_KEY) ) then
         array = 0.0_RP
         array = GetRealArrayFromString( controlVariables % StringValueForKey(GRAVITY_DIRECTION_KEY,&
                                                                             KEYWORD_LENGTH))
         if ( norm2(array(1:3)) .gt. epsilon(1.0_RP) ) then
            dimensionless_ % gravity_dir = array(1:3) / norm2(array(1:3))
         else
            dimensionless_ % gravity_dir = 0.0_RP
         end if
      end if

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
!     ------------------
!     Velocity direction
!     ------------------
!
      if ( controlVariables % ContainsKey(VELOCITY_DIRECTION_KEY) ) then
         array = 0.0_RP
         array = GetRealArrayFromString( controlVariables % StringValueForKey(VELOCITY_DIRECTION_KEY,&
                                                                             KEYWORD_LENGTH))
         dimensionless_ % vel_dir = array(1:3) / norm2(array(1:3))
      else
         print*, "*** ERROR: Introduce the 'velocity direction = [x,y,z]'"
         stop
      end if
!
!     --------------------------
!     Artificial compressibility
!     --------------------------
!
      dimensionless_ % Ma2 = POW2(refValues_ % V)/(maxval(dimensionless_ % rho) * thermodynamics_ % c02)
      dimensionless_ % invMa2 = maxval(dimensionless_ % rho) * thermodynamics_ % c02 / POW2(refValues_ % V)
!
!     ----------------
!     Density limiters
!     ----------------
!
      if ( controlVariables % ContainsKey(MAXIMUM_DENSITY_KEY) ) then
         dimensionless_ % rho_max = controlVariables % DoublePrecisionValueForKey(MAXIMUM_DENSITY_KEY) / refValues_ % rho
      else
!
!      - Default to maximum density between fluids 1 and 2
         dimensionless_ % rho_max = maxval(dimensionless_ % rho)

      end if

      if ( controlVariables % ContainsKey(MINIMUM_DENSITY_KEY) ) then
         dimensionless_ % rho_min = controlVariables % DoublePrecisionValueForKey(MINIMUM_DENSITY_KEY) / refValues_ % rho
      else
!
!      - Default to minimum density between fluids 1 and 2
         dimensionless_ % rho_min = minval(dimensionless_ % rho)

      end if
!
!     **********************************************************************
!     Set the global (proteted) thermodynamics, dimensionless, and refValues
!     **********************************************************************
!
      call setThermodynamics(thermodynamics_)
      call setDimensionless (dimensionless_ )
      call setRefValues     (refValues_     )

      END SUBROUTINE ConstructPhysicsStorage_MU
!
!     ///////////////////////////////////////////////////////
!
!     -------------------------------------------------
!!    Destructor: Does nothing for this storage
!     -------------------------------------------------
!
      SUBROUTINE DestructPhysicsStorage_MU

      END SUBROUTINE DestructPhysicsStorage_MU
!
!     //////////////////////////////////////////////////////
!
!     -----------------------------------------
!!    Descriptor: Shows the gathered data
!     -----------------------------------------
!
      SUBROUTINE DescribePhysicsStorage_MU()
         USE Headers
         use MPI_Process_Info
         IMPLICIT NONE
         real(kind=RP)  :: pRef

         if ( .not. MPI_Process % isRoot ) return

         write(STD_OUT,'(/,/)')
         call Section_Header("Loading incompressible Navier-Stokes physics")

         write(STD_OUT,'(/)')
         call SubSection_Header("Fluid data")
         write(STD_OUT,'(30X,A,A22,F10.3,A)') "->" , "Fluid 1 density: " , thermodynamics % rho(1), " kg/m^3"
         write(STD_OUT,'(30X,A,A22,F10.3,A)') "->" , "Fluid 2 density: " , thermodynamics % rho(2), " kg/m^3"
         write(STD_OUT,'(30X,A,A22,1pG10.3,A)') "->" , "Fluid 1 viscosity: " , thermodynamics % mu(1), " Pa.s"
         write(STD_OUT,'(30X,A,A22,1pG10.3,A)') "->" , "Fluid 2 viscosity: " , thermodynamics % mu(2), " Pa.s"
         write(STD_OUT,'(30X,A,A22,F10.3,A)') "->" , "Artificial compressibility c02: " , thermodynamics % c02, "(m/s)^2"

         write(STD_OUT,'(/)')
         call SubSection_Header("Reference quantities")
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference pressure: " , refValues % p, " Pa."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference density: " , refValues % rho , " kg/m^3."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference velocity: " , refValues % V , " m/s."
         write(STD_OUT,'(30X,A,A30,1pG10.3,A)') "->" , "Reference viscosity: ",refValues % mu , " Pa·s."
         write(STD_OUT,'(30X,A,A30,1pG10.3,A)') "->" , "Gravity acceleration: ",refValues % g0 , " m/s^2."

         write(STD_OUT,'(/)')
         call SubSection_Header("Dimensionless quantities")
         write(STD_OUT,'(30X,A,A20,A,F4.1,A,F4.1,A,F4.1,A)') "->" , "Velocity direction: ","[", &
                                                   dimensionless % vel_dir(1), ", ", &
                                                   dimensionless % vel_dir(2), ", ", &
                                                   dimensionless % vel_dir(3), "]"
         write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Fluid 1 Reynolds number: " , dimensionless % Re(1)
         write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Fluid 2 Reynolds number: " , dimensionless % Re(2)
         write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Froude number: " , dimensionless % Fr
         write(STD_OUT,'(30X,A,A20,A,F4.1,A,F4.1,A,F4.1,A)') "->" , "Gravity direction: ","[", &
                                                   dimensionless % gravity_dir(1), ", ", &
                                                   dimensionless % gravity_dir(2), ", ", &
                                                   dimensionless % gravity_dir(3), "]"
         write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "ACM Mach number: " , sqrt(dimensionless % Ma2)
         write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "ACM Factor: "      , dimensionless % invMa2

      END SUBROUTINE DescribePhysicsStorage_MU
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE CheckPhysics_MUInputIntegrity( controlVariables, success )
!
!        *******************************************************************
!           In this solver there are not compulsory keywords, but they
!           are given default values
!        *******************************************************************
!
         USE FTValueDictionaryClass
         USE Physics_MUKeywordsModule
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

         if ( .not. controlVariables % ContainsKey(REFERENCE_VELOCITY_KEY) ) then
            call controlVariables % AddValueForKey("1.0", REFERENCE_VELOCITY_KEY)
         end if

         if ( .not. controlVariables % ContainsKey(ARTIFICIAL_COMPRESSIBILITY_KEY) ) then
            call controlVariables % AddValueForKey("1000.0", ARTIFICIAL_COMPRESSIBILITY_KEY)
         end if

         if ( .not. controlVariables % ContainsKey(FLUID1_DENSITY_KEY)) then
            print*, "Specify density for fluid #1 using:"
            print*, "   ",trim(FLUID1_DENSITY_KEY), " = #value"
            errorMessage(STD_OUT)
            stop
         end if

         if ( .not. controlVariables % ContainsKey(FLUID2_DENSITY_KEY)) then
            print*, "Specify density for fluid #2 using:"
            print*, "   ",trim(FLUID2_DENSITY_KEY), " = #value"
            errorMessage(STD_OUT)
            stop
         end if

         if ( .not. controlVariables % ContainsKey(FLUID1_VISCOSITY_KEY)) then
            call controlVariables % AddValueForKey("0.0", FLUID1_VISCOSITY_KEY)
         end if

         if ( .not. controlVariables % ContainsKey(FLUID2_VISCOSITY_KEY)) then
            call controlVariables % AddValueForKey("0.0", FLUID2_VISCOSITY_KEY)
         end if

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
               stop

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
               stop

            else
               call controlVariables % AddValueForKey("0.0d0", GRAVITY_ACCELERATION_KEY)

            end if

         end if


      END SUBROUTINE CheckPhysics_MUInputIntegrity
!
!    **********
     END MODULE PhysicsStorage_MU
!    **********

