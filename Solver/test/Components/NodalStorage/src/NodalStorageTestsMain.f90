!
!////////////////////////////////////////////////////////////////////////
!
      PROGRAM NodalStorageTestsMain 
         USE TestSuiteManagerClass
         USE SMConstants
         IMPLICIT NONE  
         
         TYPE(TestSuiteManager) :: testSuite
         INTEGER                :: numberOfFailures
         EXTERNAL               :: testNodalStorage
         
         CALL testSuite % init()
         
         CALL testSuite % addTestSubroutineWithName(testNodalStorage,"Nodal Storage Test")
         
         CALL testSuite % performTests(numberOfFailures)
         
         CALL testSuite % finalize()
         
         IF(numberOfFailures > 0)   STOP 99

      END PROGRAM NodalStorageTestsMain
