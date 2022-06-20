!
!////////////////////////////////////////////////////////////////////////
!
      PROGRAM FacepatchTestsMain
         USE TestSuiteManagerClass
         IMPLICIT NONE 
         
         TYPE(TestSuiteManager) :: testSuite
         INTEGER                :: numberOfFailures
         
         EXTERNAL :: quadraticFaceTest
         EXTERNAL :: flatFaceTest
         EXTERNAL :: sphericalFaceTest
         EXTERNAL :: cubicFaceTest
         
         CALL testSuite % init()
         
         CALL testSuite % addTestSubroutineWithName(flatFaceTest,"Flat facePatch Tests")
         CALL testSuite % addTestSubroutineWithName(quadraticFaceTest,"Quadratic facePatch Tests")
         CALL testSuite % addTestSubroutineWithName(cubicFaceTest,"Cubic facePatch Tests")
         CALL testSuite % addTestSubroutineWithName(sphericalFaceTest,"Spherical facePatch Tests")
         
         CALL testSuite % performTests(numberOfFailures)
         CALL testSuite % finalize()
         
         IF(numberOfFailures > 0)   STOP 99
         
      END PROGRAM FacepatchTestsMain
