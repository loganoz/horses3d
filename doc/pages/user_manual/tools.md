title: Tools
---
[TOC]

# Horses tools

<!-- Advanced users can have additional control over a simulation without having to modify the source code and recompile the code. To do that, the user can provide a set of routines that are called in different stages of the simulation via the Problem file (*ProblemFile.f90*). A description of the routines of the Problem File can be found in the following section. -->

It is possible to exclusively **compile** the tools running, from the `Solver` folder:

```
make addons
```

## Stats mesh interpolation

Interpolate stats from one mesh to another. There is no need that meshes have the same geometry nor node type. The routine is completely general.

The algorithm creates a probe for each node of the output mesh `M_out`. This probe is searched in the input mesh `M_in`. Then, a high-order interpolation is computed.

**This tool should be run in sequential mode. It does not work with MPI nor OpenMP.**

### Keywords

The keywords located in the control file regarding the **input** mesh (the mesh where the fields are defined) are:

|  Keyword  |  Description  |  Value  |
|---------|-------------|---------------|
| `tool type` | Required to run this tool | `"stats mesh interpolation"`
| `stats solver` | Specifies which solver generated the stats | ns/ins/mu
| `stats geo mesh file name` | Input mesh geometry | `"path/to/mesh.msh"`
| `stats mesh file name` | Input mesh data in `hmesh` format | `"path/to/mesh.hmesh"`
| `stats qbase file name` | Stats file computed in the input mesh | `"path/to/file.stats.hsol"`
| `stats sound velocity squared file name` | In the case of Navier-Stokes simulations, stats file of the sound velocity computed in the input mesh | `"path/to/file.SoundVelocitySquared.stats.hsol"`
| `stats gradient sound velocity squared file name` | In the case of Navier-Stokes simulations, stats file of the gradient of the sound velocity computed in the input mesh | `"path/to/file.GradientSoundVelocitySquared.stats.hsol"`
| `stats lamb vector file name` | Stats file of the Lamb vector computed in the input mesh | `"path/to/file.Lamb.stats.hsol"`
| `lamb vector file name` | File of the Lamb vector computed in the input mesh | `"path/to/file.Lamb.hsol"`


The keywords located in the control file regarding the **output** mesh (the mesh where the fields should be interpolated into) are:

|  Keyword  |  Description  |  Value  |
|---------|-------------|---------------|
| `mesh file name` | Output mesh geometry | `"path/to/mesh.msh"`
| `discretization nodes` | Type of nodes of the output mesh | Gauss/Gauss-Lobatto
| `Polynomial order` | Polynomial order of the output mesh | 1

### Execution

To execute the `stats mesh interpolation` tool, it is required to set the `tool type = "stats mesh interpolation"` in the control file, and then

```
./horses.tools controlFile.control
```