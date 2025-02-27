title: Implicit Solvers with Newton Linearization
---

[TOC]

## General Keywords

The keywords for the implicit solvers are listed in the table below:

| Keyword           | Description                                                                                                                                                                                                                   | Default value |
|-------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------|
| time integration | *CHARACTER*: This is the main keyword for activating the implicit solvers. The value of it should be set to 'implicit' for the BDF solvers and to 'rosenbrock' for Rosenbrock schemes. <br> 'explicit'                                            | 'explicit'    |
| linear solver    | *CHARACTER*: Specifies the linear solver that has to be used. Options are:<br> - 'petsc': PETSc library Krylov-Subspace methods. Available in serial, but use with care (PETSc is not thread-safe, so OpenMP is not recommended). Only available in parallel (MPI) for preallocated Jacobians (see next section).<br> - 'pardiso': Intel MKL PARDISO. Only available in serial or with OpenMP.<br> - 'matrix-free gmres': A matrix-free version of the GMRES algorithm. Can be used without preconditioner or with a recursive GMRES preconditioner using 'preconditioner=GMRES'. Available in serial and parallel (OpenMP+MPI)<br> - 'smooth': Traditional iterative methods. One can select either 'smoother=WeightedJacobi' or 'smoother=BlockJacobi'.<br> - 'matrix-free smooth': A matrix-free version of the previous solver. Only available with 'smoother=BlockJacobi'. | 'petsc'       |


## Keywords for the BDF Methods

The BDF methods implemented in `HORSES3D` use a Newton's method

| Keyword                | Description                                                                                                                                                                                                                                               | Default value |
|------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------|
| bdf order              | *INTEGER*: If present, the solver uses a BDF solver of the specified order. BDF1 - BDF5 are available, and BDF2 - BDF5 require constant time steps.                                                                                                     | 1             |
| jacobian by convergence | *LOGICAL*: When .TRUE., the Jacobian is only computed when the convergence falls beneath a threshold (hard-coded). This improves performance.                                                                                                           | .FALSE.       |
| compute jacobian every | *INTEGER*: Forces the Jacobian to be computed in an interval of iterations that is specified.                                                                                                                                                           | Inf           |
| print newton info      | *LOGICAL*: If .TRUE., the information of the Newton iterations will be displayed.                                                                                                                                                                        | '.FALSE.'     |
| implicit adaptive dt   | *LOGICAL*: Specifies if the time-step should be computed according to the convergence behavior of the Newton iterative method and the linear solver.                                                                                                    | .FALSE.       |
| newton tolerance       | *REAL*: Specifies the tolerance for the Newton's method.                                                                                                                                                                                                 | \(10^{-6}\) for time-accurate simulations, or \(MaxResidual \times a\) for steady-state simulations, where \(a\) is the keyword *newton factor* |
| newton max iter        | *INTEGER*: Maximum number of Newton iterations for BDF solver.                                                                                                                                                                                          | 30            |
| linsolver max iter     | *INTEGER*: Maximum number of iterations to be taken by the linear solver. This keyword only affects iterative linear solvers.                                                                                                                            | 500           |
| newton factor          | *REAL*: In simulations that are not time-accurate, the tolerance of the Newton's method is a function of the residual: \(MaxResidual \times a\), where \(a\) is the specified value.                                                                      | \(10^{-3}\)   |
| linsolver tol factor   | *REAL*: The linear solver tolerance is a function of the absolute error of the Newton's method:  \(tol=\| e_{\infty} \|*a^i\), where \(e\) is the absolute error of the Newton's method, \(i\) is the Newton iteration number, and \(a\) is the specified value. | \(0.5\)         |
| newton first norm      | *REAL*: Specifies an assumed infinity norm of the absolute error of the Newton's method at the iteration \(0\) of the time step \(1\). This can change the behavior of the first Newton iterative method because of the dependency of the linear system tolerance on the absolute error of the Newton's method (see keyword *linsolver tol factor*).                  | \(0.2\)         |

## Keywords for the Resenbrock-Type Implicit Runge-Kutta Methods

| Keyword            | Description                                                                                                      | Default value |
|--------------------|------------------------------------------------------------------------------------------------------------------|---------------|
| rosenbrock scheme | *CHARACTER*: Rosenbrock scheme to be used. Currently, only the *RO6-6* is implemented.                         | --            |


## Jacobian Specifications
The Jacobian must be defined using a block of the form:

```markdown
#define Jacobian
   type = 2
   print info = .TRUE.
   preallocate = .TRUE.
#end
```

| Keyword      | Description                                                                                                                                                                                                                                             | Default value         |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------|
| type         | *INTEGER*: Specifies the type of Jacobian matrix to be computed. Options are:<br> 1. Numerical Jacobian: Uses a coloring algorithm and a finite difference method to compute the DG Jacobian matrix (only available with shared memory parallelization).<br> 2. Analytical Jacobian: Available with shared (OpenMP) or distributed (MPI) memory parallelization for advective and/or diffusive nonlinear conservation laws, **BUT** only for the standard DGSEM (no split-form). | **Mandatory Keyword** |
| print info   | *LOGICAL*: Specifies the verbosity of the Jacobian subroutines                                                                                                                                                                                         | .TRUE.                |
| preallocate  | *LOGICAL*: Specifies if the Jacobian must be allocated in preprocessing (.TRUE. - only available for advective/diffusive nonlinear conservation laws) or every time it is computed (.FALSE.)                                                                                                               | .FALSE.               |
