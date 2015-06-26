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
!      UserDefinedSetUp
!      UserDefinedInitialCondition(sem)
!      UserDefinedPeriodicOperation(sem)
!      UserDefinedFinalize(sem)
!      UserDefinedTermination
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedSetup  
!
!        --------------------------------
!        Called before any other routines
!        --------------------------------
!
         IMPLICIT NONE  
      END SUBROUTINE UserDefinedSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedFinalSetup(sem)
!
!        ----------------------------------------------------------------------
!        Called after the mesh is read in to allow mesh related initializations
!        or memory allocations.
!        ----------------------------------------------------------------------
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
!        ------------------------------------------------
!        Called to set the initial condition for the flow
!        ------------------------------------------------
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
!           -------------------------------------------------
!           Perturb mean flow in the expectation that it will
!           relax back to the mean flow
!           -------------------------------------------------
!
            sem % mesh % elements(eID) % Q(3,3,3,1) = 1.05_RP*sem % mesh % elements(eID) % Q(3,3,3,1)
            
         END DO 
         
      END SUBROUTINE UserDefinedInitialCondition
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedPeriodicOperation(sem, time)
!
!        ----------------------------------------------------------
!        Called at the output interval to allow periodic operations
!        to be performed
!        ----------------------------------------------------------
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
!
!        --------------------------------------------------------
!        Called after the solution computed to allow, for example
!        error tests to be performed
!        --------------------------------------------------------
!
         USE DGSEMClass
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(DGSem)  :: sem
         REAL(KIND=RP) :: time
!
!        ---------------
!        Local variables
!        ---------------
!
!         REAL(KIND=RP), ALLOCATABLE, DIMENSION(:,:,:,:) :: QExpected
!         INTEGER                                        :: N
         
         
         
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
