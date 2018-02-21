!////////////////////////////////////////////////////////////////////////
!
!      LinearSolverClass.f90
!      Created: 2017-04-12 00:20:00 +0100 
!      By: Andrés Rueda
!
!      Module to load all linear solvers
!
!////////////////////////////////////////////////////////////////////////
MODULE LinearSolverClass
   USE PetscSolverClass
   USE MKLPardisoSolverClass
   USE IterativeSolverClass
   USE MultigridSolverClass
   use MatrixFreeSmootherClass
   IMPLICIT NONE
END MODULE LinearSolverClass
