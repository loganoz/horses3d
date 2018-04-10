!
!//////////////////////////////////////////////////////
!
!   @File:    PhysicsStorage.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Wed Dec  6 17:42:24 2017
!   @Last revision date: Tue Apr 10 00:35:11 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: f29151019ab7b61620e51c8f9aaa7bca7762a0ef
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
!         character(len=KEYWORD_LENGTH), parameter    :: MOBILITY_KEY         = "mobility"
         character(len=KEYWORD_LENGTH), parameter    :: PECLET_NUMBER_KEY    = "peclet number"
         character(len=KEYWORD_LENGTH), parameter    :: INTERFACE_WIDTH_KEY  = "interface width (dimensionless)"
!         character(len=KEYWORD_LENGTH), parameter    :: INTERFACE_ENERGY_KEY = "interface energy (dimensionless)"
         CHARACTER(LEN=KEYWORD_LENGTH), DIMENSION(2) :: physicsKeywords = [INTERFACE_WIDTH_KEY, &
                                                                           PECLET_NUMBER_KEY ]
!
!        ******************
!        Optional arguments
!        ******************
!
         character(len=KEYWORD_LENGTH), parameter  :: REFERENCE_LENGTH_KEY    = "reference length"
         character(len=KEYWORD_LENGTH), parameter  :: WALL_CONTACT_ANGLE_KEY  = "wall contact angle"
!         character(len=KEYWORD_LENGTH), parameter  :: ALPHA_CONCENTRATION_KEY = "alpha concentration"
!         character(len=KEYWORD_LENGTH), parameter  :: BETA_CONCENTRATION_KEY  = "beta concentration"
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
     public    Pe
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
!    -----------------
!    The Peclet number
!    -----------------
!
     real(kind=RP), protected    :: Pe
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
!     *****************************
!     Read dimensionless properties
!     *****************************
!
      dimensionless_ % w     = controlVariables % DoublePrecisionValueForKey(INTERFACE_WIDTH_KEY)
      dimensionless_ % eps   = dimensionless_ % w 
      Pe                     = controlVariables % DoublePrecisionValueForKey(PECLET_NUMBER_KEY)
!
!     **************************************
!     Read the wall contact angle if present
!     **************************************
!
      if ( controlVariables % containsKey(WALL_CONTACT_ANGLE_KEY) ) then
         thermodynamics_ % thetaw = controlVariables % DoublePrecisionValueForKey(WALL_CONTACT_ANGLE_KEY)

      else
         thermodynamics_ % thetaw = 0.0_RP

      end if
!
!     **********************************
!     Compute the rest of the quantities
!     **********************************
!
      thermodynamics_ % rhoS  = 1.0_RP
      thermodynamics_ % M     = refValues_ % L / (thermodynamics_ % rhoS * Pe)
      thermodynamics_ % kappa = POW2(dimensionless % eps*refValues_ % L) * thermodynamics_ % rhoS
      thermodynamics_ % c_alpha = -1.0_RP
      thermodynamics_ % c_beta  =  1.0_RP

      refValues_ % time = 1.0_RP

      dimensionless_ % sigma = sqrt(2.0_RP * thermodynamics_ % kappa * thermodynamics_ % rhoS)/3.0_RP
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

         call SubSection_Header("Chemical properties")
         write(STD_OUT,'(30X,A,A40,ES10.3,A)') "->" , "Mobility: " , thermodynamics % M
         write(STD_OUT,'(30X,A,A40,ES10.3,A)') "->" , "Double-well potential height: " , thermodynamics % rhoS
         write(STD_OUT,'(30X,A,A40,ES10.3,A)') "->" , "Gradient energy coefficient: " , thermodynamics % kappa
         write(STD_OUT,'(30X,A,A40,ES10.3)') "->" , "Alpha equilibrium concentration: " , thermodynamics % c_alpha 
         write(STD_OUT,'(30X,A,A40,ES10.3)') "->" , "Beta  equilibrium concentration: " , thermodynamics % c_beta
         write(STD_OUT,'(30X,A,A40,ES10.3)') "->" , "Wall contact angle: " , thermodynamics % thetaw

      
         write(STD_OUT,'(/)')
         call SubSection_Header("Dimensionless quantities")
         write(STD_OUT,'(30X,A,A40,ES10.3)') "->" , "Interface width (dimensionless): " , dimensionless % w
         write(STD_OUT,'(30X,A,A40,ES10.3)') "->" , "Interface energy (dimensionless): " , dimensionless % sigma
         write(STD_OUT,'(30X,A,A40,ES10.3)') "->" , "Epsilon: " , dimensionless % eps

         write(STD_OUT,'(/)')
         call SubSection_Header("Reference quantities")
         write(STD_OUT,'(30X,A,A30,ES10.3,A)') "->" , "Reference length: " , refValues % L , " m."
         write(STD_OUT,'(30X,A,A30,ES10.3,A)') "->" , "Reference time: ", refValues % time, " s."
         
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

