!
!////////////////////////////////////////////////////////////////////////
!
!      FacePatchTestsMain.f90
!      Created: May 19, 2015 at 10:56 AM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      PROGRAM FacepatchTestsMain
         USE TestSuiteManagerClass
         IMPLICIT NONE 
         TYPE(TestSuiteManager) :: testSuite
         
         EXTERNAL :: quadraticFaceTest
         EXTERNAL :: flatFaceTest
         EXTERNAL :: sphericalFaceTest
         
         CALL testSuite % init()
         
         CALL testSuite % addTestSubroutineWithName(flatFaceTest,"Flat facePatch Tests")
         CALL testSuite % addTestSubroutineWithName(quadraticFaceTest,"Quadratic facePatch Tests")
         CALL testSuite % addTestSubroutineWithName(sphericalFaceTest,"Spherical facePatch Tests")
         
        CALL testSuite % performTests()
        CALL testSuite % finalize()
     
      END PROGRAM FacepatchTestsMain
