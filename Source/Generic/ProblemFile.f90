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
      IMPLICIT NONE  
      END SUBROUTINE UserDefinedSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedInitialCondition(sem)
         USE DGSEMClass
         IMPLICIT NONE
         CLASS(DGSem) :: sem
      END SUBROUTINE UserDefinedInitialCondition
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedPeriodicOperation(sem)  
         USE DGSEMClass
         IMPLICIT NONE
         CLASS(DGSem) :: sem
      END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedFinalize(sem)  
         USE DGSEMClass
         IMPLICIT NONE
         CLASS(DGSem) :: sem
      END SUBROUTINE UserDefinedFinalize
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedTermination  
         IMPLICIT NONE  
      END SUBROUTINE UserDefinedTermination
