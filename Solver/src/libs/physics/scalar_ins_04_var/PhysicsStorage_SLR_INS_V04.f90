#include "Includes.h"
      Module Physics_SLR_INS_V04KeywordsModule
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

      END MODULE Physics_SLR_INS_V04KeywordsModule
!
!////////////////////////////////////////////////////////////////////////
!
!    ******
     MODULE PhysicsStorage_SLR_INS_V04
!    ******
!
     USE SMConstants
     use FluidData_SLR_INS_V04
     use FileReadingUtilities, only: getRealArrayFromString

     IMPLICIT NONE

     private
     public    flowIsScalar_ins_v04, NCONS, NGRAD
     public    N_INS
     public    computeGradients

     public    ConstructPhysicsStorage_SLR_INS_V04, DestructPhysicsStorage_SLR_INS_V04, DescribePhysicsStorage_SLR_INS_V04
     public    CheckPhysics_SLR_INS_V04InputIntegrity
!
!    ----------------------------
!    Either NavierStokes or Euler
!    ----------------------------
!
     logical, protected :: flowIsScalar_ins_v04   = .true.
     logical, protected :: computeGradients   = .true.


!
!    --------------------------
!!   The sizes of the system
!    --------------------------
!
     INTEGER, PARAMETER :: NCONS = 7, NGRAD = 7
     INTEGER, PARAMETER :: N_INS = 3          
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
      SUBROUTINE ConstructPhysicsStorage_SLR_INS_V04( controlVariables, Lref, timeref, success )
      USE FTValueDictionaryClass
      USE Physics_SLR_INS_V04KeywordsModule
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
      CALL CheckPhysics_SLR_INS_V04InputIntegrity(controlVariables,success)
      IF(.NOT. success) RETURN
!
!     **************************************
!     Check if state gradients are requested
!     **************************************
!
      computeGradients = .true.

      END SUBROUTINE ConstructPhysicsStorage_SLR_INS_V04
!
!     ///////////////////////////////////////////////////////
!
!     -------------------------------------------------
!!    Destructor: Does nothing for this storage
!     -------------------------------------------------
!
      SUBROUTINE DestructPhysicsStorage_SLR_INS_V04

      END SUBROUTINE DestructPhysicsStorage_SLR_INS_V04
!
!     //////////////////////////////////////////////////////
!
!     -----------------------------------------
!!    Descriptor: Shows the gathered data
!     -----------------------------------------
!
      SUBROUTINE DescribePhysicsStorage_SLR_INS_V04()
         USE Headers
         use MPI_Process_Info
         IMPLICIT NONE
         real(kind=RP)  :: pRef


      END SUBROUTINE DescribePhysicsStorage_SLR_INS_V04
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE CheckPhysics_SLR_INS_V04InputIntegrity( controlVariables, success )
!
!        *******************************************************************
!           In this solver there are not compulsory keywords, but they
!           are given default values
!        *******************************************************************
!
         USE FTValueDictionaryClass
         USE Physics_SLR_INS_V04KeywordsModule
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
         



      END SUBROUTINE CheckPhysics_SLR_INS_V04InputIntegrity
!
!    **********
     END MODULE PhysicsStorage_SLR_INS_V04
!    **********
