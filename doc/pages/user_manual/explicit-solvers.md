title: Explicit Solvers
---


Explicit time integration schemes available in `HORSES3D`.
The main keywords to use it are shown in the table below

| Keyword             | Description                                                                                                                                                                                                                                                                                           | Default value |
|---------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------|
| time integration    | *CHARACTER*: This is the main keyword to activate the multigrid solvers. The value of it should be set to 'FAS' for the Full Approximation Scheme (FAS) nonlinear multigrid  solvers and to 'AnisFAS' for anisotropic FAS schemes.                                                                 | 'explicit'    |
| simulation type     | *CHARACTER*: Specifies if HORSES3D must perform a ’steady-state’ or a ’time-accurate’. If 'time-accurate' the solver switches to BDF integration and uses FAS as a pseudo problem solver. Compatible only with 'FAS'.                                                                                  | 'steady-state' |
| explicit method     | *CHARACTER*: Select desired Runge-Kutta solver. Options are: 'Euler', 'RK3', 'RK5', 'RKOpt', 'SSPRK33', 'SSPRK43, and 'Mixed RK'                                                                                                                                                                            | RK3           |
| rk order            | *INTEGER*: Order of Runge-Kutta method optimized for steady-state solver ('RKOpt'). Possible orders are from 2 to 7.                                                                                                                                                                                  | 2             |
| limit timestep      | *LOGICAL*: Activate the positivity limiter of Zhang and Shu (only for SSPRK methods).                                                                                                                                                                                                                | .false.       |
| limiter minimum     | *REAL*: Minimum value of density and pressure allowed by the limiter.                                                                                                                                                                                                                               | 1e-13         |
