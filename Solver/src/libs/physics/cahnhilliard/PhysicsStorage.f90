!
!//////////////////////////////////////////////////////
!
!   @File:    PhysicsStorage.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Wed Dec  6 17:42:24 2017
!   @Last revision date: Tue Jan 30 09:06:31 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 882f5bf78f7b6f83a4885a55d62f44ba5e21ad98
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
      Module PhysicsKeywordsModule
         IMPLICIT NONE 
         INTEGER, PARAMETER :: KEYWORD_LENGTH = 132
!
!        ******************
!        Required arguments
!        ******************
!
         character(len=KEYWORD_LENGTH), parameter    :: MOBILITY_KEY         = "mobility"
         character(len=KEYWORD_LENGTH), parameter    :: INTERFACE_WIDTH_KEY  = "interface width (dimensionless)"
         character(len=KEYWORD_LENGTH), parameter    :: INTERFACE_ENERGY_KEY = "interface energy (dimensionless)"
         CHARACTER(LEN=KEYWORD_LENGTH), DIMENSION(3) :: physicsKeywords = [MOBILITY_KEY, &
                                                                           INTERFACE_WIDTH_KEY, &
                                                                           INTERFACE_ENERGY_KEY]
!
!        ******************
!        Optional arguments
!        ******************
!
         character(len=KEYWORD_LENGTH), parameter  :: REFERENCE_LENGTH_KEY    = "reference length"
         character(len=KEYWORD_LENGTH), parameter  :: ALPHA_CONCENTRATION_KEY = "alpha concentration"
         character(len=KEYWORD_LENGTH), parameter  :: BETA_CONCENTRATION_KEY  = "beta concentration"
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
     public    N_EQN, N_GRAD_EQN, NDIM, IX, IY, IZ
     public    NCONS
     public    NGRAD
     public    NPRIM
     public    Thermodynamics, RefValues, Dimensionless
     public    Thermodynamics_t, RefValues_t, Dimensionless_t
     public    computeGradients
      
     public    ConstructPhysicsStorage, DestructPhysicsStorage, DescribePhysicsStorage
     public    CheckPhysicsInputIntegrity
!
!    ----------------------------
!    Either NavierStokes or Euler
!    ----------------------------
!
     logical, parameter :: computeGradients   = .true.
!
!    --------------------------
!!   The sizes of the NS system
!    --------------------------
!
     INTEGER, PARAMETER :: N_EQN = 1, N_GRAD_EQN = 1
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
     INTEGER, PARAMETER       :: NCONS = 1
!
!    ----------------------------------------
!!   The positions of the primitive variables
!    ----------------------------------------
!
     INTEGER, PARAMETER       :: NPRIM = 1
!
!    ---------------------------------------
!!   The positions of the gradient variables
!    ---------------------------------------
!
     INTEGER, PARAMETER  :: NGRAD = 1
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
      use Utilities, only: toLower, almostEqual
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
      type(Thermodynamics_t)           :: thermodynamics_
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
!     *************************
!     Read Reference quantities
!     *************************
!
      if ( controlVariables % containsKey(REFERENCE_LENGTH_KEY) ) then
         refValues_ % L = controlVariables % DoublePrecisionValueForKey(REFERENCE_LENGTH_KEY)

      else
         refValues_ % L = 1.0_RP       ! m

      end if
!
!     *******************************
!     Read thermodynamical properties
!     *******************************
!
      thermodynamics_ % M = controlVariables % DoublePrecisionValueForKey(MOBILITY_KEY)
   
      if ( controlVariables % containsKey(ALPHA_CONCENTRATION_KEY) ) then
         thermodynamics_ % c_alpha = controlVariables % DoublePrecisionValueForKey(ALPHA_CONCENTRATION_KEY)

      else
         thermodynamics_ % c_alpha = 0.0_RP

      end if

      if ( controlVariables % containsKey(BETA_CONCENTRATION_KEY) ) then
         thermodynamics_ % c_beta = controlVariables % DoublePrecisionValueForKey(BETA_CONCENTRATION_KEY)

      else
         thermodynamics_ % c_beta = 1.0_RP

      end if
!
!     *****************************
!     Read dimensionless properties
!     *****************************
!
      dimensionless_ % w     = controlVariables % DoublePrecisionValueForKey(INTERFACE_WIDTH_KEY)
      dimensionless_ % eps   = dimensionless_ % w / 7.071_RP
      dimensionless_ % sigma = controlVariables % DoublePrecisionValueForKey(INTERFACE_ENERGY_KEY)
!
!     **********************************
!     Compute the rest of the quantities
!     **********************************
!
      thermodynamics_ % rhoS  = (7.071_RP / 0.01508_RP ) * dimensionless_ % sigma / dimensionless_ % w
      thermodynamics_ % kappa = dimensionless_ % w * dimensionless_ % sigma * POW2(refValues_ % L) / (7.071_RP * 0.01508_RP)

      refValues_ % time = POW2(refValues_ % L) / thermodynamics_ % M
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

         if ( .not. MPI_Process % isRoot ) return 

         call Section_Header("Loading Cahn-Hilliard physics")

         write(STD_OUT,'(/,/)')

         call SubSection_Header("Reference quantities")
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reynolds length: " , refValues % L , " m."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference time: ", refValues % time, " s."
         
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

