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
!//////////////////////////////////////////////////////////////////////// 
! 
      MODULE UserDefinedFunctions
      
      CONTAINS 
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
            USE FTValueDictionaryClass
            IMPLICIT NONE
            CLASS(DGSem)            :: sem
            TYPE(FTValueDictionary) :: controlVariables
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedInitialCondition(sem , controlVariables)
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
            
            TYPE(DGSem)              :: sem
            class(FTValueDictionary) :: controlVariables
            EXTERNAL                 :: initialStateSubroutine
                     
            INTEGER     :: i, j, k, eID
            
            DO eID = 1, SIZE(sem % mesh % elements)
               DO k = 0, sem % mesh % elements(eID) % Nxyz(3)
                  DO j = 0, sem % mesh % elements(eID) % Nxyz(2)
                     DO i = 0, sem % mesh % elements(eID) % Nxyz(1)
                        CALL UniformFlowState( sem % mesh % elements(eID) % geom % x(:,i,j,k), 0.0_RP, &
                                               sem % mesh % elements(eID) % Q(i,j,k,1:N_EQN) )
                                                     
                     END DO
                  END DO
               END DO 
!
!              -------------------------------------------------
!              Perturb mean flow in the expectation that it will
!              relax back to the mean flow
!              -------------------------------------------------
!
               sem % mesh % elements(eID) % Q(1,1,1,1) = 1.05_RP*sem % mesh % elements(eID) % Q(1,1,1,1)
               
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
            USE BoundaryConditionFunctions
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
!~            INTEGER                            :: numberOfFailures
!~            CHARACTER(LEN=29)                  :: testName           = "27 element uniform flow tests"
!~            REAL(KIND=RP)                      :: maxError
!~            REAL(KIND=RP), ALLOCATABLE         :: QExpected(:,:,:,:)
!~            INTEGER                            :: eID
!~            INTEGER                            :: i, j, k, N
!~            TYPE(FTAssertionsManager), POINTER :: sharedManager
!~!
!~!           -----------------------------------------------------------------------
!~!           Expected Values. Note they will change if the run parameters change and
!~!           when the eigenvalue computation for the time step is fixed. These 
!~!           results are for the Mach 0.5 and rusanov solvers.
!~!           -----------------------------------------------------------------------
!~!
!~            INTEGER                            :: expectedIterations      = 45
!~            REAL(KIND=RP)                      :: expectedResidual        = 7.8469500886160862E-011
            
!~            CALL initializeSharedAssertionsManager
!~            sharedManager => sharedAssertionsManager()
            
!~            N = sem % mesh % elements(1) % Nxyz(1) ! This works here because all the elements have the same order
!~            CALL FTAssertEqual(expectedValue = expectedIterations, &
!~                               actualValue   =  sem % numberOfTimeSteps, &
!~                               msg           = "Number of time steps to tolerance")
!~            CALL FTAssertEqual(expectedValue = expectedResidual, &
!~                               actualValue   = sem % maxResidual, &
!~                               tol           = 1.d-3, &
!~                               msg           = "Final maximum residual")
            
!~            ALLOCATE(QExpected(0:N,0:N,0:N,N_EQN))
            
!~            maxError = 0.0_RP
!~            DO eID = 1, SIZE(sem % mesh % elements)
!~               DO k = 0, sem % mesh % elements(eID) % Nxyz(3)
!~                  DO j = 0, sem % mesh % elements(eID) % Nxyz(2)
!~                     DO i = 0, sem % mesh % elements(eID) % Nxyz(1)
!~                        CALL UniformFlowState( sem % mesh % elements(eID) % geom % x(:,i,j,k), 0.0_RP, &
!~                                               QExpected(i,j,k,1:N_EQN) )
!~                     END DO
!~                  END DO
!~               END DO
!~               maxError = MAXVAL(ABS(QExpected - sem % mesh % elements(eID) % Q))
!~            END DO
!~            CALL FTAssertEqual(expectedValue = 0.0_RP, &
!~                               actualValue   = maxError, &
!~                               tol           = 1.d-10, &
!~                               msg           = "Maximum error")
            
            
!~            CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)
   
!~            IF ( sharedManager % numberOfAssertionFailures() == 0 )     THEN
!~               WRITE(6,*) testName, " ... Passed"
!~            ELSE
!~               WRITE(6,*) testName, " ... Failed"
!~               WRITE(6,*) "NOTE: Failure is expected when the max eigenvalue procedure is fixed."
!~               WRITE(6,*) "      When that is done, re-compute the expected values and modify this procedure"
!~               error stop 99
!~            END IF 
!~            WRITE(6,*)
            
!~            CALL finalizeSharedAssertionsManager
!~            CALL detachSharedAssertionsManager
            
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
      
!
!=====================================================================================================
!=====================================================================================================
!
!