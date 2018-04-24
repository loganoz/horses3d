!
!//////////////////////////////////////////////////////
!
!   @File:    PhysicsStorage_CH.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Thu Apr 19 17:24:30 2018
!   @Last revision date: Tue Apr 24 11:13:20 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 60804273321199c0675663c7d4f1c517987552a7
!
!//////////////////////////////////////////////////////
!
!
!//////////////////////////////////////////////////////
!
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
      Module Physics_CHKeywordsModule
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
!         character(len=KEYWORD_LENGTH), parameter    :: INTERFACE_ENERGY_KEY = "interface energy (multiphase)"
         CHARACTER(LEN=KEYWORD_LENGTH), DIMENSION(2) :: physics_CHKeywords = [INTERFACE_WIDTH_KEY, &
                                                                              PECLET_NUMBER_KEY ]
!
!        ******************
!        Optional arguments
!        ******************
!
         character(len=KEYWORD_LENGTH), parameter  :: WALL_CONTACT_ANGLE_KEY  = "wall contact angle"
!         character(len=KEYWORD_LENGTH), parameter  :: ALPHA_CONCENTRATION_KEY = "alpha concentration"
!         character(len=KEYWORD_LENGTH), parameter  :: BETA_CONCENTRATION_KEY  = "beta concentration"
      END MODULE Physics_CHKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
!    
!    ******
     MODULE PhysicsStorage_CH
!    ******
!
     USE SMConstants
     USE Physics_CHKeywordsModule
     use FluidData_CH
     
     IMPLICIT NONE

     private
     public    NCOMP

     public    ConstructPhysicsStorage_CH, DestructPhysicsStorage_CH, DescribePhysicsStorage_CH
     public    CheckPhysicsCHInputIntegrity

     integer, parameter    :: NCOMP = 1
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
      SUBROUTINE ConstructPhysicsStorage_CH( controlVariables, Lref, tref, success )
      USE FTValueDictionaryClass
      use Utilities, only: toLower, almostEqual
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(FTValueDictionary) :: controlVariables
      real(kind=RP),    intent(inout)    :: Lref
      real(kind=RP),    intent(inout)    :: tref
      LOGICAL                 :: success
!
!     ---------------
!     Local variables
!     ---------------
!
      CHARACTER(LEN=KEYWORD_LENGTH) :: keyword
      type(Multiphase_t)            :: multiphase_
!
!     --------------------
!     Collect input values
!     --------------------
!
      success = .TRUE.
      CALL CheckPhysicsCHInputIntegrity(controlVariables,success)
      IF(.NOT. success) RETURN 
!
!     *****************************
!     Read multiphase properties
!     *****************************
!
      multiphase_ % w   = controlVariables % DoublePrecisionValueForKey(INTERFACE_WIDTH_KEY)
      multiphase_ % eps = multiphase_ % w
      multiphase_ % Pe  = controlVariables % DoublePrecisionValueForKey(PECLET_NUMBER_KEY)
!
!     **************************************
!     Read the wall contact angle if present
!     **************************************
!
      if ( controlVariables % containsKey(WALL_CONTACT_ANGLE_KEY) ) then
         multiphase_ % thetaw = controlVariables % DoublePrecisionValueForKey(WALL_CONTACT_ANGLE_KEY)

      else
         multiphase_ % thetaw = 0.0_RP

      end if
!
!     **********************************
!     Compute the rest of the quantities
!     **********************************
!
      multiphase_ % rhoS    = 1.0_RP
      multiphase_ % M       = Lref / (multiphase_ % rhoS * multiphase_ % Pe)
      multiphase_ % kappa   = POW2(multiphase_ % eps*Lref) * multiphase_ % rhoS
      multiphase_ % c_alpha = -1.0_RP
      multiphase_ % c_beta  = 1.0_RP

      multiphase_ % sigma   = sqrt(2.0_RP * multiphase_ % kappa * multiphase_ % rhoS)/3.0_RP
!
!     ************************************
!     Set the global (proteted) multiphase
!     ************************************
!
      call setMultiphase( multiphase_ )

      CALL DescribePhysicsStorage_CH()

      END SUBROUTINE ConstructPhysicsStorage_CH
!
!     ///////////////////////////////////////////////////////
!
!     -------------------------------------------------
!!    Destructor: Does nothing for this storage
!     -------------------------------------------------
!
      SUBROUTINE DestructPhysicsStorage_CH
      
      END SUBROUTINE DestructPhysicsStorage_CH
!
!     //////////////////////////////////////////////////////
!
!     -----------------------------------------
!!    Descriptor: Shows the gathered data
!     -----------------------------------------
!
      SUBROUTINE DescribePhysicsStorage_CH()
         USE Headers
         use MPI_Process_Info
         IMPLICIT NONE

         if ( .not. MPI_Process % isRoot ) return 

         call Section_Header("Loading Cahn-Hilliard physics")

         write(STD_OUT,'(/,/)')

         call SubSection_Header("Chemical properties")
         write(STD_OUT,'(30X,A,A40,ES10.3,A)') "->" , "Mobility: " , multiphase % M
         write(STD_OUT,'(30X,A,A40,ES10.3,A)') "->" , "Double-well potential height: " , multiphase % rhoS
         write(STD_OUT,'(30X,A,A40,ES10.3,A)') "->" , "Gradient energy coefficient: " , multiphase % kappa
         write(STD_OUT,'(30X,A,A40,ES10.3)') "->" , "Alpha equilibrium concentration: " , multiphase % c_alpha 
         write(STD_OUT,'(30X,A,A40,ES10.3)') "->" , "Beta  equilibrium concentration: " , multiphase % c_beta
         write(STD_OUT,'(30X,A,A40,ES10.3)') "->" , "Wall contact angle: " , multiphase % thetaw

      
         write(STD_OUT,'(/)')
         call SubSection_Header("Dimensionless quantities")
         write(STD_OUT,'(30X,A,A40,ES10.3)') "->" , "Interface width (dimensionless): " , multiphase % w
         write(STD_OUT,'(30X,A,A40,ES10.3)') "->" , "Interface energy (dimensionless): " , multiphase % sigma
         write(STD_OUT,'(30X,A,A40,ES10.3)') "->" , "Epsilon: " , multiphase % eps

      END SUBROUTINE DescribePhysicsStorage_CH
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE CheckPhysicsCHInputIntegrity( controlVariables, success )  
         USE FTValueDictionaryClass
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
         
         DO i = 1, SIZE(physics_CHKeywords)
            obj => controlVariables % objectForKey(physics_CHKeywords(i))
            IF ( .NOT. ASSOCIATED(obj) )     THEN
               PRINT *, "Input file is missing entry for keyword: ",physics_CHKeywords(i)
               success = .FALSE. 
            END IF  
         END DO  
         
      END SUBROUTINE CheckPhysicsCHInputIntegrity
!
!    **********       
     END MODULE PhysicsStorage_CH
!    **********

