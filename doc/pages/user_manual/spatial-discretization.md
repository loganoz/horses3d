title: Spatial Discretization
---

`HORSES3D` uses a high-order discontinuous Galerkin Spectral Element (DGSEM) method. As part of the numerical method, several options of volume and surface discretization options can be selected. Some options are only available in certain solvers.

| Keyword | Description | Default value |
|---------|-------------|---------------|
| Discretization nodes | *CHARACTER*: Node type to use in the solution. Options are: <br> - Gauss <br> - Gauss-Lobatto | 'Gauss' |
| Inviscid discretization | *CHARACTER*: DG discretization of inviscid fluxes. Options are: <br> - Standard <br> - Split-Form | 'Standard' |
| Averaging | *CHARACTER*: Type of averaging function to use in numerical fluxes and split-forms (if in use). Options are: <br> - Standard <br> - Morinishi <br> - Ducros <br> - Kennedy-Gruber <br> - Pirozzoli <br> - Entropy conserving <br> - Chandrasekar <br> - Skew-symmetric 1 (only for Incompressible Navier–Stokes) <br> - Skew-symmetric 2 (only for Incompressible Navier–Stokes) | -- |
| Viscous discretization | *CHARACTER*: Method for viscous fluxes. Options are: <br> - BR1 <br> - BR2 <br> - IP | 'BR2' |
| Riemann solver | *CHARACTER*: Riemann solver for inviscid fluxes. Options are: <br> - Central <br> - Roe (Only for compressible Navier–Stokes) <br> - Standard Roe (Only for compressible Navier–Stokes) <br> - Roe-Pike (Only for compressible Navier–Stokes) <br> - Low dissipation Roe (Only for compressible Navier–Stokes) <br> - Lax-Friedrichs (Only for compressible and Incompressible Navier–Stokes) <br> - ES Lax-Friedrichs (Only for compressible Navier–Stokes, not including RANS) <br> - u-diss (Only for compressible Navier–Stokes, not including RANS) <br> - Rusanov (Only for compressible Navier–Stokes) <br> - Matrix dissipation (Only for compressible Navier–Stokes) <br> - Viscous NS (Only for compressible Navier–Stokes, not including RANS) <br> - Exact (Only for Incompressible Navier-Stokes, and Multiphase) | **Mandatory keyword** |



