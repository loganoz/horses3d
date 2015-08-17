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
         SUBROUTINE UserDefinedFinalSetup(sem)
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            USE DGSEMClass
            IMPLICIT NONE
            CLASS(DGSem) :: sem
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedInitialCondition(sem)
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
            
            TYPE(DGSem)      :: sem
            EXTERNAL         :: initialStateSubroutine
                     
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
               sem % mesh % elements(eID) % Q(3,3,3,1) = 1.05_RP*sem % mesh % elements(eID) % Q(3,3,3,1)
               
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

            write(6,*) "This test case has no expected solution yet"
            
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
      SUBROUTINE ExternalGradientForBoundaryName( x, t, nHat, U_x, U_y, U_z, boundaryType )
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
      REAL(KIND=RP)   , INTENT(INOUT) :: U_x(N_GRAD_EQN), U_y(N_GRAD_EQN), U_z(N_GRAD_EQN)
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryType
!
!     ---------------
!     Local variables
!     ---------------
!
      IF ( boundarytype == "freeslipwall" )                   THEN
         CALL FreeSlipNeumann( x, t, nHat, U_x, U_y, U_z )
      ELSE IF ( boundaryType == "noslipadiabaticwall" )       THEN 
         CALL  NoSlipAdiabaticWallNeumann( x, t, nHat, U_x, U_y, U_z)
      ELSE IF ( boundarytype == "noslipisothermalwall")       THEN 
         CALL NoSlipIsothermalWallNeumann( x, t, nHat, U_x, U_y, U_z )
      ELSE
         CALL UniformFlowNeumann( x, t, nHat, U_x, U_y, U_z )
      END IF

      END SUBROUTINE ExternalGradientForBoundaryName

