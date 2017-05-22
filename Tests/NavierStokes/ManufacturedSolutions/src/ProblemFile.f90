!
!////////////////////////////////////////////////////////////////////////
!
!      ProblemFile.f90
!      Created: June 26, 2015 at 8:47 AM 
!      By: David Kopriva  
!
!      The Problem File contains user defined procedures
!      that are used to "personalize" i.e. define a specific
!      problem to be solved. These procedures include initial conditions,
!      exact solutions (e.g. for tests), etc. and allow modifications 
!      without having to modify the main code.
!
!      The procedures, *even if empty* that must be defined are
!
!      UserDefinedStartup
!      UserDefinedInitialCondition(sem)
!      UserDefinedPeriodicOperation(sem)
!      UserDefinedFinalize(sem)
!      UserDefinedTermination
!
!      *** This problem file sets up a subsonic point source *** 
!
!//////////////////////////////////////////////////////////////////////// 
!
      MODULE UserDefinedDataStorage
         USE SMConstants
         IMPLICIT NONE 
         REAL(KIND=RP) :: rad0, f, h 
      END MODULE UserDefinedDataStorage
!
!//////////////////////////////////////////////////////////////////////// 
! 
MODULE UserDefinedFunctions
   USE SMConstants
   IMPLICIT NONE
   
   CHARACTER(LEN=LINE_LENGTH) :: ManSolType
!
!     ========      
      CONTAINS
!     ========
!
         SUBROUTINE UserDefinedStartup  
!
!        --------------------------------
!        Called before any other routines
!        --------------------------------
!
            IMPLICIT NONE  
         END SUBROUTINE UserDefinedStartup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalSetup(sem, controlVariables)
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            USE DGSEMClass
            IMPLICIT NONE
            CLASS(DGSem)             :: sem
            class(FTValueDictionary) :: controlVariables
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedInitialCondition(sem, controlVariables)
!
!           ------------------------------------------------
!           Called to set the initial condition for the flow
!           ------------------------------------------------
!
            USE SMConstants
            USE DGSEMClass
            USE PhysicsStorage
            USE BoundaryConditionFunctions
            IMPLICIT NONE
            
            TYPE(DGSem)                :: sem
            class(FTValueDictionary)   :: controlVariables
            !----------------------------------------------
            ABSTRACT INTERFACE
               SUBROUTINE GaussianPerturbSub(x,Q)
                  USE PhysicsStorage
                  USE SMConstants
                  REAL(KIND=RP) :: x(3)
                  REAL(KIND=RP) :: Q(N_EQN)
               END SUBROUTINE GaussianPerturbSub
            END INTERFACE
            PROCEDURE(GaussianPerturbSub), POINTER :: GaussianPerturb
            INTEGER                    :: i, j, k, eID
            
            
            ManSolType = controlVariables % StringValueForKey("manufactured solution",LINE_LENGTH)
            SELECT CASE (ManSolType)
               CASE('3D')
                  GaussianPerturb => GaussianPerturbUnitCube
               CASE('2D')
                  GaussianPerturb => GaussianPerturbUnitSquare
               CASE DEFAULT
                  print*, 'ERROR: Not recognized manufactured solution type "', TRIM(ManSolType), '"'
                  STOP
            END SELECT
            
            DO eID = 1, SIZE(sem % mesh % elements)
               DO k = 0, sem % mesh % elements(eID) % Nxyz(3)
                  DO j = 0, sem % mesh % elements(eID) % Nxyz(2)
                     DO i = 0, sem % mesh % elements(eID) % Nxyz(1)
                        CALL ManufacturedSolutionState( sem % mesh % elements(eID) % geom % x(:,i,j,k), 0.0_RP, &
                                               sem % mesh % elements(eID) % Q(i,j,k,1:N_EQN) )    !ZeroFlowState !ManufacturedSolutionState !UniformFlowState
                        
                        CALL GaussianPerturb (sem % mesh % elements(eID) % geom % x(:,i,j,k),  &
                                              sem % mesh % elements(eID) % Q(i,j,k,1:N_EQN)) 
                     END DO
                  END DO
               END DO 
!
!              -------------------------------------------------
!              Perturb mean flow in the expectation that it will
!              relax back to the mean flow
!              -------------------------------------------------
!
!~                sem % mesh % elements(eID) % Q(3,3,1,1) = 1.05_RP*sem % mesh % elements(eID) % Q(3,3,1,1)
               
            END DO 
            
         END SUBROUTINE UserDefinedInitialCondition
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedPeriodicOperation(sem, time)
!
!           ----------------------------------------------------------
!           Called at the output interval to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            USE DGSEMClass
            IMPLICIT NONE
            CLASS(DGSem)  :: sem
            REAL(KIND=RP) :: time
            
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(sem, time)
            USE Physics
            USE FTAssertions
!
!           --------------------------------------------------------
!           Called after the solution computed to allow, for example
!           error tests to be performed
!           --------------------------------------------------------
!
            USE DGSEMClass
            IMPLICIT NONE
!
!           ---------
!           Arguments
!           ---------
!
            CLASS(DGSem)  :: sem
            REAL(KIND=RP) :: time
!
!           ---------------
!           Local variables
!           ---------------
!
            CHARACTER(LEN=29)                  :: testName           = "Box Around Cyrcle test"
            REAL(KIND=RP)                      :: maxError
            REAL(KIND=RP), ALLOCATABLE         :: QExpected(:,:,:,:)
            INTEGER                            :: eID
            INTEGER                            :: i, j, k
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
!
!           -----------------------------------------------------------------------------------------
!           Expected solutions. 
!           front 0.0 manufacturedsol
!           back 0.0 manufacturedsol
!           bottom 0.0 manufacturedsol   / freeslipWall
!           top 0.0 manufacturedsol      / freeslipWall
!           left 0.0 manufacturedsol
!           right 0.0 manufacturedsol
!           -----------------------------------------------------------------------------------------
!
!
!           ------------------------------------------------
!           Expected Solutions:
!           Number of iterations are for CFL of 0.4 (0.3 for NS 3D), for
!           the rusanov solver and mach = 1.5, N = 6 in the needed directions (3D/2D)
!           ------------------------------------------------
!
            INTEGER        :: iterations
            REAL(KIND=RP)  :: residuals
!
            IF (flowIsNavierStokes) THEN
               SELECT CASE (ManSolType)
                  CASE('3D')
                     iterations = 5387
                     residuals  = 9.9708685752375459E-011
                  CASE('2D')
                     iterations = 7328
                     residuals  = 9.9887209614735184E-011
               END SELECT
            ELSE
               SELECT CASE (ManSolType)
                  CASE('3D')
                     iterations = 3660
                     residuals  = 9.9688257648722356E-011
                  CASE('2D')
                     iterations = 12824
                     residuals  = 9.9886321436315484E-011
               END SELECT
            END IF
            
            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()
            
            CALL FTAssertEqual(expectedValue = iterations, &
                               actualValue   =  sem % numberOfTimeSteps, &
                               msg           = "Number of time steps to tolerance")
            CALL FTAssertEqual(expectedValue = residuals, &
                               actualValue   = sem % maxResidual, &
                               tol           = 1.d-3, &
                               msg           = "Final maximum residual")
            
            !ALLOCATE(QExpected(0:sem % spA % N,0:sem % spA % N,0:sem % spA % N,N_EQN))
            
            ! maxError = 0.0_RP
            ! DO eID = 1, SIZE(sem % mesh % elements)
            !    DO k = 0, sem % spA % N
            !       DO j = 0, sem % spA % N
            !          DO i = 0, sem % spA % N 
            !             CALL pointSourceFlowSolution( sem % mesh % elements(eID) % geom % x(:,i,j,k), &
            !                                           QExpected(i,j,k,1:N_EQN), success )
            !          END DO
            !       END DO
            !    END DO
            !    maxError = MAXVAL(ABS(QExpected - sem % mesh % elements(eID) % Q))
            ! END DO
            ! CALL FTAssertEqual(expectedValue = ERRORs(N), &
            !                    actualValue   = maxError, &
            !                    tol           = 1.d-5, &
            !                    msg           = "Maximum error")
            
            
            CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)
   
            IF ( sharedManager % numberOfAssertionFailures() == 0 )     THEN
               WRITE(6,*) testName, " ... Passed"
               WRITE(6,*) "This test case has no expected solution yet, only checks the residual after 100 iterations."
            ELSE
               WRITE(6,*) testName, " ... Failed"
               WRITE(6,*) "NOTE: Failure is expected when the max eigenvalue procedure is changed."
               WRITE(6,*) "      If that is done, re-compute the expected values and modify this procedure"
               STOP 99
            END IF 
            WRITE(6,*)
            
            CALL finalizeSharedAssertionsManager
            CALL detachSharedAssertionsManager
            
         END SUBROUTINE UserDefinedFinalize
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedTermination
!
!        -----------------------------------------------
!        Called at the the end of the main driver after 
!        everything else is done.
!        -----------------------------------------------
!
         IMPLICIT NONE  
      END SUBROUTINE UserDefinedTermination
      
END MODULE UserDefinedFunctions

