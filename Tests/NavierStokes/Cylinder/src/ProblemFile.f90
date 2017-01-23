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
            
            TYPE(DGSem)              :: sem
            class(FTValueDictionary) :: controlVariables
            EXTERNAL                 :: initialStateSubroutine
                     
            INTEGER     :: i, j, k, eID
            
            DO eID = 1, SIZE(sem % mesh % elements)
               DO k = 0, sem % spA % N
                  DO j = 0, sem % spA % N
                     DO i = 0, sem % spA % N 
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
!               sem % mesh % elements(eID) % Q(3,3,3,1) = 1.05_RP*sem % mesh % elements(eID) % Q(3,3,3,1)
               
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
            INTEGER                            :: i, j, k, N
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
!
!           -----------------------------------------------------------------------------------------
!           Expected solutions. 
!           InnerCylinder 0.0 NoSlipAdiabaticWall
!           Front 0.0 Inflow
!           bottom 0.0 FreeSlipWall
!           top 0.0 FreeSlipWall
!           Back 0.0 Inflow
!           Left 0.0 Inflow
!           Right 0.0 OutflowSpecifyP 
!           -----------------------------------------------------------------------------------------
!
!
!           ------------------------------------------------
!           Expected Solutions: Wall conditions on the sides
!           Number of iterations are for CFL of 0.3, for
!           the roe solver and mach = 0.3
!           ------------------------------------------------
!
            INTEGER                            :: iterations(3:7) = [100, 0, 0, 0, 0]
            REAL(KIND=RP), DIMENSION(3:7)      :: residuals = [239.99443625202642, 0E-011, &
                                                               0E-011, 0E-011, &
                                                               0E-011]
!
            N = sem % spA % N
            
            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()
            
            CALL FTAssertEqual(expectedValue = iterations(N), &
                               actualValue   =  sem % numberOfTimeSteps, &
                               msg           = "Number of time steps to tolerance")
            CALL FTAssertEqual(expectedValue = residuals(N), &
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
!
!=====================================================================================================
!=====================================================================================================
!
!
      SUBROUTINE externalStateForBoundaryName( x, t, nHat, Q, boundaryType )
!
!     ----------------------------------------------
!     Set the boundary conditions for the mesh by
!     setting the external state for each boundary.
!     ----------------------------------------------
!
      USE BoundaryConditionFunctions
      USE UserDefinedDataStorage
      USE MeshTypes
      
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: Q(N_EQN)
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryType
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP)   :: pExt
      LOGICAL         :: success

      IF ( boundarytype == "freeslipwall" )             THEN
         CALL FreeSlipWallState( x, t, nHat, Q )
      ELSE IF ( boundaryType == "noslipadiabaticwall" ) THEN 
         CALL  NoSlipAdiabaticWallState( x, t, Q)
      ELSE IF ( boundarytype == "noslipisothermalwall") THEN 
         CALL NoSlipIsothermalWallState( x, t, Q )
      ELSE IF ( boundaryType == "outflowspecifyp" )     THEN 
         pExt =  ExternalPressure()
         CALL ExternalPressureState ( x, t, nHat, Q, pExt )
      ELSE 
         CALL UniformFlowState( x, t, Q ) 
      END IF

      END SUBROUTINE externalStateForBoundaryName
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ExternalGradientForBoundaryName( x, t, nHat, GradU, boundaryType )
!
!     ------------------------------------------------
!     Set the boundary conditions for the mesh by
!     setting the external gradients on each boundary.
!     ------------------------------------------------
!
      USE BoundaryConditionFunctions
      USE MeshTypes
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: GradU(3,N_GRAD_EQN)
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryType
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP) :: U_x(N_GRAD_EQN), U_y(N_GRAD_EQN), U_z(N_GRAD_EQN)

      U_x(:) = GradU(1,:)
      U_y(:) = GradU(2,:)
      U_z(:) = GradU(3,:)

      IF ( boundarytype == "freeslipwall" )                   THEN
         CALL FreeSlipNeumann( x, t, nHat, U_x, U_y, U_z )
      ELSE IF ( boundaryType == "noslipadiabaticwall" )       THEN 
         CALL  NoSlipAdiabaticWallNeumann( x, t, nHat, U_x, U_y, U_z )
      ELSE IF ( boundarytype == "noslipisothermalwall")       THEN 
         CALL NoSlipIsothermalWallNeumann( x, t, nHat, U_x, U_y, U_z )
      ELSE
         CALL UniformFlowNeumann( x, t, nHat, U_x, U_y, U_z )
      END IF

      GradU(1,:) = U_x(:)
      GradU(2,:) = U_y(:)
      GradU(3,:) = U_z(:)

      END SUBROUTINE ExternalGradientForBoundaryName

