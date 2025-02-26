!////////////////////////////////////////////////////////////////////////
!
!      Module to load all linear solvers
!
!////////////////////////////////////////////////////////////////////////
module LinearSolverClass
#if defined(SCALAR_INS_V04)
   use MKLPardisoSolverClass
#else
   use PetscSolverClass
   use MKLPardisoSolverClass
   use StaticCondensationSolverClass
   use IterativeSolverClass
   use LinearMultigridSolverClass
   use MatrixFreeSmootherClass
   use MatrixFreeGMRESClass
#endif
   implicit none
end module LinearSolverClass