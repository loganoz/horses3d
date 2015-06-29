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
      MODULE UserDefinedFunctions
      
      CONTAINS 
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
      
