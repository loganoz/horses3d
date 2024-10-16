#include "Includes.h"
      Module Physics_SLRKeywordsModule
         IMPLICIT NONE
         INTEGER, parameter :: KEYWORD_LENGTH = 132


         !PARTICLES 
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

      END MODULE Physics_SLRKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
!
!    ******
     MODULE PhysicsStorage_SLR
!    ******
!
     USE SMConstants
     use FluidData_SLR
     use FileReadingUtilities, only: getRealArrayFromString

     IMPLICIT NONE

     private
     public    NCONS, NGRAD
     public    computeGradients

     public    ConstructPhysicsStorage_SLR, DestructPhysicsStorage_SLR, DescribePhysicsStorage_SLR
     public    CheckPhysics_SLRInputIntegrity
!
!    ----------------------------
!    Either NavierStokes or Euler
!    ----------------------------
!
     logical, protected :: computeGradients   = .true.

!
!    --------------------------
!!   The sizes of the system
!    --------------------------
!
     INTEGER, PARAMETER :: NCONS = 1, NGRAD = 1
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
      SUBROUTINE ConstructPhysicsStorage_SLR( controlVariables, Lref, timeref, success )
      USE FTValueDictionaryClass
      USE Physics_SLRKeywordsModule
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
      CALL CheckPhysics_SLRInputIntegrity(controlVariables,success)
      IF(.NOT. success) RETURN
!
!     **************************************
!     Check if state gradients are requested
!     **************************************
!
      computeGradients = .true.

      END SUBROUTINE ConstructPhysicsStorage_SLR
!
!     ///////////////////////////////////////////////////////
!
!     -------------------------------------------------
!!    Destructor: Does nothing for this storage
!     -------------------------------------------------
!
      SUBROUTINE DestructPhysicsStorage_SLR

      END SUBROUTINE DestructPhysicsStorage_SLR
!
!     //////////////////////////////////////////////////////
!
!     -----------------------------------------
!!    Descriptor: Shows the gathered data
!     -----------------------------------------
!
      SUBROUTINE DescribePhysicsStorage_SLR()
         USE Headers
         use MPI_Process_Info
         IMPLICIT NONE
         real(kind=RP)  :: pRef


      END SUBROUTINE DescribePhysicsStorage_SLR
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE CheckPhysics_SLRInputIntegrity( controlVariables, success )
!
!        *******************************************************************
!           In this solver there are not compulsory keywords, but they
!           are given default values
!        *******************************************************************
!
         USE FTValueDictionaryClass
         USE Physics_SLRKeywordsModule
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
         



      END SUBROUTINE CheckPhysics_SLRInputIntegrity
!
!    **********
     END MODULE PhysicsStorage_SLR
!    **********
