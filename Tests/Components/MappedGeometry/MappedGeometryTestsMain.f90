!
!////////////////////////////////////////////////////////////////////////
!
!      MappedGeometryTestsMain.f90
!      Created: May 22, 2015 at 12:46 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      PROGRAM MappedGeometryTestsMain 
         USE TestSuiteManagerClass
         IMPLICIT NONE
         
         TYPE(TestSuiteManager) :: testSuite
         EXTERNAL               :: cubeTest
         EXTERNAL               :: cylinderTestGeometry
        
         CALL testSuite % init()
         
         CALL testSuite % addTestSubroutineWithName(cubeTest,testName = "Cube Geometry")
         CALL testSuite % addTestSubroutineWithName(cylinderTestGeometry,testName = "Cylindrical Geometry")
         
         CALL testSuite % performTests()
         CALL testSuite % finalize()
         
         
      END PROGRAM MappedGeometryTestsMain
