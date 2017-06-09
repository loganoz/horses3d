!
!////////////////////////////////////////////////////////////////////////
!
!      NSLite3D.f90
!      Created: May 21, 2015 at 12:56 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE TestDivergence
      USE FTAssertions
      USE DGSEMClass
      USE DGTimeDerivativeMethods
      USE ProlongToFacesProcedures
      USE SetupModule
      USE DGViscousDiscretization
      
      IMPLICIT NONE
!
!     ------------
!     Declarations
!     ------------
!
      INTEGER            :: eID
      INTEGER            :: nElement
      REAL(KIND=RP)      :: maxE
      LOGICAL            :: success
      CHARACTER(LEN=132) :: msg
      CHARACTER(LEN=132), EXTERNAL :: lastPathComponent
      INTEGER            :: N(3)
!
!     ------------------------------
!     Read in the mesh for this test
!     ------------------------------
!
      CALL setUpDGSEM(meshFileName = meshFileNames(testFileCount), &
                      success = success)
      msg = "Read in mesh " // lastPathComponent(meshFileNames(testFileCount))
      CALL FTAssert(success,msg)
      IF(.NOT.success) RETURN 
!
!     -----------------
!     Perform the tests
!     -----------------
!
!
!     ----------------------------------------------------
!     In this point, do not use any viscous discretization
!     ----------------------------------------------------
!
      deallocate( ViscousMethod )
      nElement =  SIZE(sem % mesh % elements)
      DO eID = 1, nElement
         N = sem % mesh % elements(eID) % Nxyz
         CALL ProlongToFaces(sem % mesh % elements(eId), sem % spA(N(1),N(2),N(3)))
      END DO
      
      CALL computeRiemannFluxes(sem,0.0_RP)
      
      call TimeDerivative_ComputeQDot( sem % mesh , sem % spA , 0.0_RP )
!
!     ------------------------------------------------
!     Check the divergence of the different components
!     ------------------------------------------------
!
      DO eID = 1, nElement
          WRITE(msg,'(A,I3)') "Gradient of F = x on element ",eID
          maxE = MAXVAL(ABS(sem % mesh % elements(eID) % QDot(:,:,:,1)+1.0_RP))
          CALL FTAssertEqual(expectedValue = 0.0_RP, &
                             actualValue = maxE,     &
                             tol = 1.e-9_RP,            &
                             msg = msg)
                             
          WRITE(msg,'(A,I3)') "Gradient of F = y on element ",eID
          maxE = MAXVAL(ABS(sem % mesh % elements(eID) % QDot(:,:,:,2)+1.0_RP))
          CALL FTAssertEqual(expectedValue = 0.0_RP, &
                             actualValue = maxE,     &
                             tol = 1.e-9_RP,            &
                             msg = msg)
                             
          WRITE(msg,'(A,I3)') "Gradient of F = z on element ",eID
          maxE = MAXVAL(ABS(sem % mesh % elements(eID) % QDot(:,:,:,3)+1.0_RP))
          CALL FTAssertEqual(expectedValue = 0.0_RP, &
                             actualValue = maxE,     &
                             tol = 1.e-9_RP,            &
                             msg = msg)
                             
          WRITE(msg,'(A,I3)') "Gradient of F = const on element ",eID
          maxE = MAXVAL(ABS(sem % mesh % elements(eID) % QDot(:,:,:,4)))
          CALL FTAssertEqual(expectedValue = 0.0_RP, &
                             actualValue = maxE,     &
                             tol = 1.e-9_RP,            &
                             msg = msg)
                             
          WRITE(msg,'(A,I3)') "Gradient of F = x + y + z on element ",eID
          maxE = MAXVAL(ABS(sem % mesh % elements(eID) % QDot(:,:,:,5)+3.0_RP))
          CALL FTAssertEqual(expectedValue = 0.0_RP, &
                             actualValue = maxE,     &
                             tol = 1.e-9_RP,            &
                             msg = msg)
          
      END DO 
!
!     -------------------------------------------
!     Destroy the mesh in preparation for another
!     -------------------------------------------
!
      CALL sem % destruct()
      testFileCount = testFileCount + 1
      IF(testFileCount == SIZE(meshFileNames)+1) testFileCount = 1
      
      END SUBROUTINE TestDivergence      
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE TestGradients
      USE FTAssertions
      USE DGSEMClass
      USE DGTimeDerivativeMethods
      USE ProlongToFacesProcedures
      USE SetupModule
      
      IMPLICIT NONE
!
!     ------------
!     Declarations
!     ------------
!
      INTEGER            :: eID
      INTEGER            :: nElement
      REAL(KIND=RP)      :: maxE
      LOGICAL            :: success
      CHARACTER(LEN=132) :: msg
      CHARACTER(LEN=132), EXTERNAL :: lastPathComponent
      INTEGER            :: N(3)
!
!     ------------------------------
!     Read in the mesh for this test
!     ------------------------------
!
      CALL setUpDGSEM(meshFileName = meshFileNames(testFileCount), &
                      success = success)
      msg = "Read in mesh " // lastPathComponent(meshFileNames(testFileCount))
      CALL FTAssert(success,msg)
      IF(.NOT.success) RETURN 

!     -----------------
!     Perform the tests
!     -----------------
!
      nElement =  SIZE(sem % mesh % elements)
      DO eID = 1, nElement
         N = sem % mesh % elements(eID) % Nxyz
         CALL ProlongToFaces(sem % mesh % elements(eId), sem % spA(N(1),N(2),N(3)))
      END DO
      
      IF ( flowIsNavierStokes )     THEN
         CALL DGSpatial_ComputeGradient( sem % mesh , sem % spA , 0.0_RP , sem % externalState , sem % externalGradients ) 
      END IF

      CALL computeRiemannFluxes(sem,0.0_RP)
      
      call TimeDerivative_ComputeQDot( sem % mesh , sem % spA , 0.0_RP )
!
!     ------------------------------------------------
!     Check the divergence of the different components
!     ------------------------------------------------
!
      DO eID = 1, nElement
          WRITE(msg,'(A,I3)') "Gradient of F = x on element ",eID
          maxE = MAXVAL(ABS(sem % mesh % elements(eID) % QDot(:,:,:,1)+1.0_RP))
          CALL FTAssertEqual(expectedValue = 0.0_RP, &
                             actualValue = maxE,     &
                             tol = 1.e-9_RP,         &
                             msg = msg)
                             
          WRITE(msg,'(A,I3)') "Gradient of F = y on element ",eID
          maxE = MAXVAL(ABS(sem % mesh % elements(eID) % QDot(:,:,:,2)+1.0_RP))
          CALL FTAssertEqual(expectedValue = 0.0_RP, &
                             actualValue = maxE,     &
                             tol = 1.e-9_RP,         &
                             msg = msg)
                             
          WRITE(msg,'(A,I3)') "Gradient of F = z on element ",eID
          maxE = MAXVAL(ABS(sem % mesh % elements(eID) % QDot(:,:,:,3)+1.0_RP))
          CALL FTAssertEqual(expectedValue = 0.0_RP, &
                             actualValue = maxE,     &
                             tol = 1.e-9_RP,         &
                             msg = msg)
                             
          WRITE(msg,'(A,I3)') "Gradient of F = const on element ",eID
          maxE = MAXVAL(ABS(sem % mesh % elements(eID) % QDot(:,:,:,4)))
          CALL FTAssertEqual(expectedValue = 0.0_RP, &
                             actualValue = maxE,     &
                             tol = 1.e-9_RP,         &
                             msg = msg)
                             
          WRITE(msg,'(A,I3)') "Gradient of F = x + y + z on element ",eID
          maxE = MAXVAL(ABS(sem % mesh % elements(eID) % QDot(:,:,:,5)+3.0_RP))
          CALL FTAssertEqual(expectedValue = 0.0_RP, &
                             actualValue = maxE,     &
                             tol = 1.e-9_RP,         &
                             msg = msg)
          
      END DO 
!
!     -------------------------------------------
!     Destroy the mesh in preparation for another
!     -------------------------------------------
!
      CALL sem % destruct()
      testFileCount = testFileCount + 1
      IF(testFileCount == SIZE(meshFileNames)+1) testFileCount = 1
      
      END SUBROUTINE TestGradients
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE TestInterpolationToFaces
      USE setupModule
      USE ProlongToFacesProcedures
      USE FTAssertions
      USE DGSEMClass
      
      IMPLICIT NONE
!
!     ------------
!     Declarations
!     ------------
!
      LOGICAL            :: success
      CHARACTER(LEN=132) :: msg
      CHARACTER(LEN=132), EXTERNAL :: lastPathComponent
!
!     ------------------------------
!     Read in the mesh for this test
!     ------------------------------
!
      CALL setUpDGSEM(meshFileName = meshFileNames(testFileCount), &
                      success = success)
      msg = "Read in mesh " // lastPathComponent(meshFileNames(testFileCount))
      CALL FTAssert(success,msg)
      IF(.NOT.success) RETURN 
!
!     -----------------
!     Conduct the tests
!     -----------------
!
      CALL interpolateToFaces(sem)
!
!     -------------------------------------------
!     Destroy the mesh in preparation for another
!     -------------------------------------------
!
      CALL sem % destruct()
      testFileCount = testFileCount + 1
      IF(testFileCount == SIZE(meshFileNames)+1) testFileCount = 1
      
      END SUBROUTINE TestInterpolationToFaces
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SetInitialCondition( sem, initialStateSubroutine )
         USE SMConstants
         USE DGSEMClass
         USE PhysicsStorage
         IMPLICIT NONE
         
         TYPE(DGSem)      :: sem
         EXTERNAL         :: initialStateSubroutine
                  
         INTEGER     :: i, j, k, eID
         INTEGER     :: N(3)
         
         DO eID = 1, SIZE(sem % mesh % elements)
            N = sem % mesh % elements(eID) % Nxyz
            DO k = 0, N(3)
               DO j = 0, N(2)
                  DO i = 0, N(1) 
                     CALL initialStateSubroutine( sem % mesh % elements(eID) % geom % x(:,i,j,k), 0.0_RP, &
                                                  sem % mesh % elements(eID) % Q(i,j,k,1:N_EQN) )
                  END DO
               END DO
            END DO 
         END DO 
         
      END SUBROUTINE SetInitialCondition
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE interpolateToFaces(sem)
         USE FTAssertions
         USE DGSEMClass
         USE PhysicsStorage
         USE ProlongToFacesProcedures
         USE DGTimeDerivativeMethods
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(DGSem)  :: sem
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER           :: eID
         INTEGER           :: i, j
         INTEGER           :: fce
         INTEGER           :: N(3)
         REAL(KIND=RP)     :: x(3), Qexpected(N_EQN), Qactual(N_EQN), emax
         CHARACTER(LEN=72) :: msg
!
!        ---------------------------------------------------------------------
!        Interpolate to the faces. Then compare to the expected values at the 
!        boundary points. For this test the interpolations should be exact.
!        ---------------------------------------------------------------------
!
         
         DO eID = 1, SIZE(sem % mesh % elements)
            N = sem % mesh % elements(eID) % Nxyz
            CALL ProlongToFaces(sem % mesh % elements(eId), sem % spA(N(1),N(2),N(3)))
            DO fce = 1, 6
               emax = 0.0_RP
               DO j = 0, N(axisMap(2,fce))
                  DO i = 0, N(axisMap(1,fce))
!
!                    --------------
!                    Expected value
!                    --------------
!
                     x = sem % mesh % elements(eID) % geom % xb(:,i,j,fce)
                     CALL initialFlowState(x, 0.0_RP, Qexpected)
!
!                    ------------
!                    Actual value
!                    ------------
!
                     Qactual = sem % mesh % elements(eId) % Qb(:,i,j,fce)
                     emax = MAX(MAXVAL(ABS(Qactual-Qexpected)),emax)
                        
                  END DO
               END DO   
                        
               WRITE(msg,'(A,i1,A,i3)') "Face values of solution on face ", fce, " in element ", eID
               CALL FTAssertEqual(expectedValue = 0.0_RP,actualValue = emax,tol = 1.0e-11_RP,msg = msg)
            END DO
         END DO
         
      END SUBROUTINE interpolateToFaces
