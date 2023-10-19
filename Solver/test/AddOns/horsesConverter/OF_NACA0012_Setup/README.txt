This case is an example to run OpenFOAM with HORSES3D mesh and convert VTK result of OF to HORSES3D solution file (.hsol)

Note: You need OpenFOAM software 

1. This setup requires mesh files (.hmesh and .bmesh) from NACA0012 test.
2. The Horses mesh files is then converted to OF polyMesh files
3. Copy polyMesh folder to constant/ folder
4. Modify the boundary type of boundary file (patch, wall or symmetry) according to actual boundary condition.
5. To run OpenFOAM see run_OF.sh
6. Convert the result to VTK with foamToVTK (ascii format)
7. To convert VTK result to .hsol see horsesConverter.convert
8. See User Manual for more details