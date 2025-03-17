!////////////////////////////////////////////////////////////////////////
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
   implicit none
end module LinearSolverClass