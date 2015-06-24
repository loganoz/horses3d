!
!////////////////////////////////////////////////////////////////////////
!
!      Setup.f90
!      Created: June 19, 2015 at 12:52 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      MODULE setupModule
      USE ControlVariablesModule
      USE DGSEMClass
      
      TYPE( NSLiteControlVariables )  :: controlVariables
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
!
!     ---------------
!     Initializations
!     ---------------
!
      CALL constructSharedBCModule
      CALL ReadInputFile( controlVariables )
      CALL ConstructPhysicsStorage( controlVariables % mach,               &
                                    controlVariables % RE,                 &
                                    0.72_RP,                               &
                                    controlVariables % flowIsNavierStokes )
!
!     --------------------------
!     Get the file names to test
!     --------------------------
!
      fUnit = UnusedUnit()
      OPEN(UNIT = fUnit, FILE = controlVariables % inputFileName)
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
         
         EXTERNAL            :: initialFlowState
         EXTERNAL            :: externalBoundaryState, externalGradientState
!
!        ----------------
!        Set up the DGSEM
!        ----------------
!
         CALL ConstructDGSem(self              = sem, &
                             polynomialOrder   = controlVariables % polynomialOrder,&
                             meshFileName      = meshFileName,  &
                             externalState     = externalBoundaryState,  &
                             externalGradients = externalGradientState,  &
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
