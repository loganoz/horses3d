!
!////////////////////////////////////////////////////////////////////////
!
!   @File:    TestSuiteModule.f90
!   @Author:  David Kopriva
!   @Created: February 21, 2013 11:21 AM
!   @Last revision date: Tue Sep 28 11:34:35 2021
!   @Last revision author: Wojciech Laskowski (wj.laskowski@upm.es)
!   @Last revision commit: 0b8b49ef742bce3e02d3138ef5e95597b5d3a726
!
!////////////////////////////////////////////////////////////////////////
!
!> The TestSuiteManager class defines methods to easily
!> put together and run a suite of unit tests. 
!>
!>
!> The tests are managed by an instance of the
!>**TestSuiteManager** class. It is designed to be used with minimal fuss. You
!>
!>- Initialize the test suite
!>- Add test subroutines
!>- Have the testSuiteManager perform the tests
!>- Finalize the test suite manager
!>
!># Usage: #
!>
!>##Definition
!>
!>      TYPE(TestSuiteManager) :: testSuite
!>
!>##Initialization
!>         call testSuite % init()
!>
!>##Creating a test ###
!>
!>   A test is a subroutine with interface
!>
!>         ABSTRACT INTERFACE
!>            SUBROUTINE testSuiteSubroutine()
!>            END SUBROUTINE testSuiteSubroutine
!>         END INTERFACE
!>   
!>   that (typically) includes unit test calls. You add
!>   a test suite function by the add subroutine
!>   
!>         CALL testSuite % addTestSubroutineWithName(SubroutineName, description)
!>
!>   where 
!>
!> - SubroutineName = a subroutine with the interface as above, and 
!> - description = a CHARACTER(LEN=128) character string that names the test
!>   
!>##Setting the output location ###
!>   Set the unit to which the output is written by
!>
!>         CALL testSuite % setOutputUnit(iUnit)
!>
!>##Running tests ###
!>   To run the tests call
!>
!>         CALL testSuite % performTests() OR
!>         CALL testSuite % performTests(numFailed)
!>   
!>##Finalizing the test suite ###
!>   When done, call
!>
!>         CALL testSuite % finalize()
!
!////////////////////////////////////////////////////////////////////////
!
      Module TestSuiteManagerClass
      USE FTAssertions
      IMPLICIT NONE
      PRIVATE 
            
      ABSTRACT INTERFACE
         SUBROUTINE testSuiteFunction()
         END SUBROUTINE testSuiteFunction
      END INTERFACE

      TYPE TestCaseRecord
         LOGICAL                                       :: passed
         CHARACTER(LEN=128)                            :: testName
         TYPE(FTAssertionsManager)   , POINTER         :: assertionsManager
         PROCEDURE(testSuiteFunction), POINTER, NOPASS :: TestSubroutine
         TYPE(TestCaseRecord), POINTER                 :: next
      END TYPE TestCaseRecord
      
      TYPE, PUBLIC :: TestSuiteManager
         INTEGER                       :: numberOfTests
         INTEGER                       :: stdOut = 6
         TYPE(TestCaseRecord), POINTER :: testCasesHead => NULL()
         TYPE(TestCaseRecord), POINTER :: testCasesTail => NULL()
         CONTAINS
         PROCEDURE :: init     => initializeTestSuiteManager
         PROCEDURE :: finalize => finalizeTestSuiteManager
         PROCEDURE :: addTestSubroutineWithName
         PROCEDURE :: performTests
         PROCEDURE :: setOutputUnit
      END TYPE TestSuiteManager
!
!     ========      
      CONTAINS
!     ========      
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE initializeTestSuiteManager(self)
         IMPLICIT NONE
         CLASS(TestSuiteManager) :: self

            self % testCasesHead => NULL()
            self % testCasesTail => NULL()
            self % numberOfTests = 0
            
      END SUBROUTINE initializeTestSuiteManager
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE setOutputUnit(self,iUnit)
         IMPLICIT NONE 
         CLASS(TestSuiteManager) :: self
         INTEGER                 :: iUnit
         self % stdOut = iUnit
      END SUBROUTINE setOutputUnit    
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE addTestSubroutineWithName(self,testSubroutine ,testName)
         IMPLICIT NONE
         CLASS(TestSuiteManager)       :: self
         EXTERNAL                      :: testSubroutine
         CHARACTER(LEN=*)              :: testName
         TYPE(TestCaseRecord), POINTER :: newTestCase
         
         INTERFACE
            SUBROUTINE testSubroutine()
            END SUBROUTINE  testSubroutine
         END INTERFACE
         
         ALLOCATE(newTestCase)
         newTestCase % testName     = TRIM(ADJUSTL(testName))
         newTestCase % TestSubroutine => testSubroutine
         newTestCase % next         => NULL()
         newTestCase % passed       = .TRUE.
         self % numberOfTests       = self % numberOfTests + 1
         
         IF ( ASSOCIATED(self % testCasesHead) )     THEN
            self % testCasesTail % next => newTestCase
            self % testCasesTail      => newTestCase 
         ELSE
            self % testCasesHead => newTestCase
            self % testCasesTail => newTestCase
         END IF 
         
      END SUBROUTINE addTestSubroutineWithName    
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE finalizeTestSuiteManager(self)
         IMPLICIT NONE
         CLASS(TestSuiteManager)       :: self
         TYPE(TestCaseRecord), POINTER :: tmp, current
         
         IF ( .NOT.ASSOCIATED(self % testCasesHead) )     THEN
           RETURN 
         END IF 
         
         current => self % testCasesHead
         DO WHILE (ASSOCIATED(current))
            tmp => current % next
            
            IF(ASSOCIATED(current % assertionsManager)) THEN
               DEALLOCATE(current % assertionsManager)
            END IF 
            
            DEALLOCATE(current)
            current => tmp
         END DO

         self % testCasesHead => NULL()
         self % testCasesTail => NULL()
         self % numberOfTests = 0
         
      END SUBROUTINE finalizeTestSuiteManager
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE performTests(self, numberOfFailedTestsRet)  
          IMPLICIT NONE  
!
!         ---------
!         Arguments
!         ---------
!
          CLASS(TestSuiteManager)             :: self
          INTEGER                  , OPTIONAL :: numberOfFailedTestsRet
!
!         ---------------
!         Local variables
!         ---------------
!
          TYPE(TestCaseRecord)     , POINTER  :: current
          TYPE(FTAssertionsManager), POINTER  :: sharedManager
          INTEGER                             :: maxMessageLength, numberOfFailedTests
          
          numberOfFailedTests = 0
          maxMessageLength    = 0
          
          WRITE(self % stdOut,*)
          WRITE(self % stdOut,*) "                   ////////////////////////////////"
          WRITE(self % stdOut,*) "                   ////    Begin Test Suites   ////"
          WRITE(self % stdOut,*) "                   ////////////////////////////////"
          WRITE(self % stdOut,*)
          
          current => self % testCasesHead
          DO WHILE (ASSOCIATED(current))
          
            CALL initializeSharedAssertionsManager
            sharedManager               => sharedAssertionsManager()
            current % assertionsManager => sharedManager
            
            CALL current % TestSubroutine
            
            IF ( sharedManager % numberOfAssertionFailures() /= 0 )     THEN
               numberOfFailedTests = numberOfFailedTests + 1 
               current % passed = .FALSE.
            END IF 
               
            CALL sharedManager % SummarizeAssertions(current % testName,self % stdOut)
            CALL detachSharedAssertionsManager
            
            maxMessageLength = MAX(maxMessageLength,LEN_TRIM(current % testName))
            
            current => current % next
          END DO
        
          WRITE(self % stdOut,*)
          WRITE(self % stdOut,*) "   **********************************************************"
          WRITE(self % stdOut,*) "                     Summary of failed test suites:"
          WRITE(self % stdOut,'(i6,A,i3)')  numberOfFailedTests," suite(s) failed out of ", self % numberOfTests 
          WRITE(self % stdOut,*) "   **********************************************************"
          
          WRITE(self % stdOut,*)
          WRITE(self % stdOut,*) "                   ////////////////////////////////////"
          WRITE(self % stdOut,*) "                   ////    Test Suites Completed   ////"
          WRITE(self % stdOut,*) "                   ////////////////////////////////////"
          WRITE(self % stdOut,*)
        
!
!         ------------------
!         Test matrix output
!         ------------------
!
          WRITE(self % stdOut,*)
          WRITE(self % stdOut,*) "////////////////////////////////"
          WRITE(self % stdOut,*) "////   Test Status Matrix   ////"
          WRITE(self % stdOut,*) "////////////////////////////////"
          WRITE(self % stdOut,*)
        
          current => self % testCasesHead
          DO WHILE (ASSOCIATED(current))
            
            IF ( current % passed )     THEN
               WRITE(self % stdOut,*) current % testName(1:maxMessageLength), " ... Passed"
            ELSE 
               WRITE(self % stdOut,*) current % testName(1:maxMessageLength), " ... F A I L E D"
            END IF 
            
            current => current % next
          END DO
        
          IF(PRESENT(numberOfFailedTestsRet)) numberOfFailedTestsRet = numberOfFailedTests
          
      END SUBROUTINE performTests    
      
      END Module TestSuiteManagerClass    