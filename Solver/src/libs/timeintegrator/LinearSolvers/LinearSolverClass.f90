!////////////////////////////////////////////////////////////////////////
!
!      LinearSolverClass.f90
!      Created: 2017-04-12 00:20:00 +0100 
!      By: Andrés Rueda
!
!      Module to load all linear solvers
!
!////////////////////////////////////////////////////////////////////////
module LinearSolverClass
   use PetscSolverClass
   use MKLPardisoSolverClass
   use StaticCondensationSolverClass
   use IterativeSolverClass
   use LinearMultigridSolverClass
   use MatrixFreeSmootherClass
   use MatrixFreeGMRESClass
   use MAGMASolverClass
   use CuSparseSolverClass
   implicit none
end module LinearSolverClass
