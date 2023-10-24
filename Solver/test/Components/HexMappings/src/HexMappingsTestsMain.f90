!
!////////////////////////////////////////////////////////////////////////
!
!      HexMappingsTests.f90
!      Created: May 19, 2015 at 4:37 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      PROGRAM HexMappingsTests
         USE TestSuiteManagerClass
         IMPLICIT NONE
         TYPE(TestSuiteManager) :: testSuite
         INTEGER                :: numberOFFailures
         EXTERNAL               :: testCube
         EXTERNAL               :: testGenHexAsCube
         EXTERNAL               :: testFlaredHex8
         EXTERNAL               :: testFlaredHex8AsGeneral
         EXTERNAL               :: testCubeWithHexMapper 
         EXTERNAL               :: testFlaredHex8AsGeneralWithMapper
         EXTERNAL               :: testCylindricalMappingWithMapper
      
         CALL testSuite % init()
         
         CALL testSuite % addTestSubroutineWithName(testCube,"Cube Hex8")
         CALL testSuite % addTestSubroutineWithName(testGenHexAsCube,"Cube as a general hex")
         CALL testSuite % addTestSubroutineWithName(testFlaredHex8,"Flared hex as Hex8")
         CALL testSuite % addTestSubroutineWithName(testFlaredHex8AsGeneral,"Flared hex as General Hex8")
         CALL testSuite % addTestSubroutineWithName(testCubeWithHexMapper,"Cube Hex8 using mapper")
         CALL testSuite % addTestSubroutineWithName(testFlaredHex8AsGeneralWithMapper, "Flared hex as General Hex8 with mapper")
         CALL testSuite % addTestSubroutineWithName(testCylindricalMappingWithMapper,"Cylindrical Geometry with mapper")
        
         CALL testSuite % performTests(numberOfFailures)
         
         CALL testSuite % finalize()
         
         IF(numberOFFailures > 0) ERROR STOP 99
         
      
      END PROGRAM HexMappingsTests