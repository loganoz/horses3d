!
!////////////////////////////////////////////////////////////////////////
!
      MODULE setupModule
      USE FTValueDictionaryClass
      USE DGSEMClass
      USE mainKeywordsModule
      
      TYPE( FTValueDictionary )       :: controlVariables
      TYPE(DGsem)                     :: sem
      INTEGER                         :: testFileCount = 1
      CHARACTER(LEN=132), ALLOCATABLE :: meshFileNames(:)
      
      CONTAINS 
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE setup  
      USE SMConstants
      USE PhysicsStorage
      USE SharedBCModule
      use Utilities, only: UnusedUnit
      
      IMPLICIT NONE
      INTEGER           :: fUnit
      INTEGER           :: numberOfMeshfiles
      INTEGER           :: j
      REAL(KIND=RP)     :: machVal, reVal
!
!     ---------------
!     Initializations
!     ---------------
!
      CALL constructSharedBCModule
      CALL controlVariables % initWithSize(16)
      CALL ReadInputFile( controlVariables )
      machVal = controlVariables % quadValueForKey(machNumberKey)
      reVal   = controlVariables % quadValueForKey(reynoldsNumberKey)
      
      CALL ConstructPhysicsStorage( machVal,     &
                                    reVal,       &
                                    0.72_RP,     &
                                    controlVariables % logicalValueForKey(flowIsNavierStokesKey) )
!
!     --------------------------
!     Get the file names to test
!     --------------------------
!
      fUnit = UnusedUnit()
      OPEN(UNIT = fUnit, FILE = controlVariables % stringValueForKey(meshFileNameKey,KEYWORD_LENGTH))
         READ(fUnit,*) numberOfMeshfiles
         ALLOCATE(meshFileNames(numberOfMeshfiles))
         DO j = 1, numberOfMeshfiles
            READ(fUnit,'(A)') meshFileNames(j)
         END DO  
      CLOSE(fUnit)
      testFileCount = 1
      
      END SUBROUTINE setup
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE setUpDGSEM(meshFileName, success)  
         IMPLICIT NONE  
!
!        ---------
!        Arguments
!        ---------
!
         CHARACTER(LEN=132)  :: meshFileName
         LOGICAL             :: success
         INTEGER             :: N(3)
         
         EXTERNAL            :: initialFlowState
         EXTERNAL            :: externalBoundaryState, externalGradientState
!
!        ----------------
!        Set up the DGSEM
!        ----------------
!
         N = controlVariables % integerValueForKey(polynomialOrderKey)
         call InitializeNodalStorage(maxval(N))
         
         CALL ConstructDGSem(self              = sem                  , &
                             polynomialOrder   = N                    , &
                             controlVariables  = controlVariables     , &
                             meshFileName_     = meshFileName         , &
                             externalState     = externalBoundaryState, &
                             externalGradients = externalGradientState, &
                             success           = success)
         IF(.NOT. success)     RETURN 
!
!        ----------------------
!        Set the initial values
!        ----------------------
!
         CALL SetInitialCondition( sem, initialFlowState )
         
      END SUBROUTINE setUpDGSEM
      
      END MODULE setupModule
