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
   use ConjugateGradientClass
   implicit none
end module LinearSolverClass