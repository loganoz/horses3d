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
      use NodalStorageClass
      use FileReaders      , only: ReadControlFile
      use Utilities        , only: UnusedUnit
      
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
      CALL ReadControlFile( controlVariables )
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
      use NodalStorageClass
      use SpatialDiscretization
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
         type(BCFunctions_t)                 :: BCFunctions(1)
         procedure(BCState_FCN)              :: externalBoundaryState
         procedure(BCGradients_FCN)          :: ExternalGradientState
!
!        ----------------
!        Set up the DGSEM
!        ----------------
!
         N = controlVariables % integerValueForKey(polynomialOrderKey)
         call InitializeNodalStorage(controlVariables, maxval(N))
         
         BCFunctions(1) % externalState     => externalBoundaryState
         BCFunctions(1) % externalGradients => externalGradientState
         CALL ConstructDGSem(self              = sem                  , &
                             polynomialOrder   = N                    , &
                             controlVariables  = controlVariables     , &
                             meshFileName_     = meshFileName         , &
                             BCFunctions       = BCFunctions          , &
                             success           = success)
         IF(.NOT. success)     RETURN 
!
!        ----------------------
!        Set the initial values
!        ----------------------
!
         CALL SetInitialCondition( sem, initialFlowState )
   
         call Initialize_SpaceAndTimeMethods(controlVariables, sem % mesh)
         
      END SUBROUTINE setUpDGSEM
      
      END MODULE setupModule
