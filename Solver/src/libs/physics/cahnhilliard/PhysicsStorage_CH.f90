!
!//////////////////////////////////////////////////////
!
!   @File:    PhysicsStorage_CH.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Thu Apr 19 17:24:30 2018
!   @Last revision date: Thu Jul 26 17:26:21 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: ba557cd23630b1bd1f528599b9b33812f58d1f7b
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
         character(len=KEYWORD_LENGTH), parameter    :: MOBILITY_KEY     = "mobility"
         character(len=KEYWORD_LENGTH), parameter    :: INTERFACE_WIDTH_KEY   = "interface width"
         character(len=KEYWORD_LENGTH), parameter    :: INTERFACE_TENSION_KEY = "interface tension"
         CHARACTER(LEN=KEYWORD_LENGTH), DIMENSION(3) :: physics_CHKeywords = [INTERFACE_WIDTH_KEY, &
                                                                              MOBILITY_KEY,   &
                                                                              INTERFACE_TENSION_KEY ]
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

      multiphase_ = Multiphase_t()
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
      multiphase_ % eps   = controlVariables % DoublePrecisionValueForKey(INTERFACE_WIDTH_KEY)
      multiphase_ % M     = controlVariables % DoublePrecisionValueForKey(MOBILITY_KEY)
      multiphase_ % sigma = controlVariables % DoublePrecisionValueForKey(INTERFACE_TENSION_KEY)
!
!     ************************************
!     Set the global (proteted) multiphase
!     ************************************
!
      call setMultiphase( multiphase_ )

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

      
         write(STD_OUT,'(/)')
         call SubSection_Header("Dimensionless quantities")
         write(STD_OUT,'(30X,A,A40,ES10.3)') "->" , "Interface tension: " , multiphase % sigma
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

