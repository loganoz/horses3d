!
!////////////////////////////////////////////////////////////////////////
!
!      TestDivergenceMain.f90
!      Created: June 18, 2015 at 12:44 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      PROGRAM TestDivergenceMain
      USE FTAssertions
      USE TestSuiteManagerClass
      USE setupModule
      IMPLICIT NONE
      
      TYPE(TestSuiteManager) :: testSuite
      
      EXTERNAL                     :: TestDivergence
      EXTERNAL                     :: TestInterpolationToFaces
      CHARACTER(LEN=132), EXTERNAL :: lastPathComponent
      CHARACTER(LEN=132)           :: msg
      INTEGER                      :: j
      
      CALL testSuite % init()
!
!     ---------------------
!     Set up the mesh, etc.
!     ---------------------
!
      CALL setup
!
!     -----------------------------------------------------------------------
!     For good measure, test that the interpolations to the faces are correct
!     -----------------------------------------------------------------------
!
      DO j = 1, SIZE(meshFileNames)
         msg = "Test face interpolations on mesh " // lastPathComponent(meshFileNames(j)) 
         CALL testSuite % addTestSubroutineWithName(testSubroutine = TestInterpolationToFaces, &
                                                    testName       = msg)
      END DO
!
!     -----------------------------------------
!     Test the divergence on each of the meshes
!     -----------------------------------------
!
      DO j = 1, SIZE(meshFileNames)
         msg = "Test divergence computation on mesh " // lastPathComponent(meshFileNames(j)) 
         CALL testSuite % addTestSubroutineWithName(testSubroutine = TestDivergence, &
                                                    testName       = msg)
      END DO
      
      CALL testSuite % performTests()
      CALL testSuite % finalize()
      
      END PROGRAM TestDivergenceMain
!
!//////////////////////////////////////////////////////////////////////// 
! 
      FUNCTION lastPathComponent(path)  
         IMPLICIT NONE  
         CHARACTER(LEN=132) :: path
         CHARACTER(LEN=132) :: lastPathComponent
         INTEGER            :: l, strLen
         
         strLen = LEN_TRIM(path)
         l      = INDEX(STRING = path, SUBSTRING = "/", BACK = .TRUE.)
         lastPathComponent = path(l+1:strLen)
          
      END FUNCTION lastPathComponent
