!
!////////////////////////////////////////////////////////////////////////
!
!      Setup.f90
!      Created: June 19, 2015 at 12:52 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      Module mainKeywordsModule
         IMPLICIT NONE 
         INTEGER, PARAMETER :: KEYWORD_LENGTH = 132
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: machNumberKey           = "mach number"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: reynoldsNumberKey       = "reynolds number"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: aoaThetaKey             = "aoa theta"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: aoaPhiKey               = "aoa phi"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: flowIsNavierStokesKey   = "flowisnavierstokes"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: polynomialOrderKey      = "polynomial order"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: cflKey                  = "cfl"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: meshFileNameKey         = "mesh file name"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: restartKey              = "restart"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: restartFileNameKey      = "restart file name"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: numberOfTimeStepsKey    = "number of time steps"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: outputIntervalKey       = "output interval"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: convergenceToleranceKey = "convergence tolerance"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: numberOfPlotPointsKey   = "number of plot points"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: numberOfBoundariesKey   = "number of boundaries"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: plotFileNameKey         = "plot file name"
         CHARACTER(LEN=KEYWORD_LENGTH), DIMENSION(15) :: mainKeywords =  [machNumberKey,           &
                                                                          reynoldsNumberKey,       &
                                                                          aoaThetaKey,             &
                                                                          aoaPhiKey,               &
                                                                          flowIsNavierStokesKey,   &
                                                                          polynomialOrderKey,      &
                                                                          cflKey,                  &
                                                                          meshFileNameKey,         &
                                                                          restartKey,              &
                                                                          restartFileNameKey,      &
                                                                          numberOfTimeStepsKey,    &
                                                                          outputIntervalKey,       &
                                                                          convergenceToleranceKey, &
                                                                          numberOfPlotPointsKey,   &
                                                                          plotFileNameKey]
      END MODULE mainKeywordsModule
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
      
      IMPLICIT NONE
      INTEGER, EXTERNAL :: UnusedUnit
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
         CALL ConstructDGSem(self              = sem                  , &
                             polynomialOrder   = N                    , &
                             meshFileName      = meshFileName         , &
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
