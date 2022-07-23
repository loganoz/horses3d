## HORSES3D High-Order (DG) Spectral Element Solver                  

MIT License                                      
Copyright (c) 2021 NUMATH https://numath.dmae.upm.es                 

## Synopsis
-----------
**HORSES3D** is a multiphysics environment where the compressible Navier-Stokes equations, the incompressible Navier–Stokes equations, the Cahn–Hilliard equation and entropy–stable variants are solved. Arbitrary high–order, p–anisotropic discretisations are used, including static and dynamic p–adaptation methods (feature-based and truncation error-based). Explicit and implicit time-steppers for steady and time-marching solutions are available, including efficient multigrid and preconditioners. Numerical and analytical Jacobian computations with a coloring algorithm have been implemented. Multiphase flows are solved using a diffuse interface model: Navier–Stokes/Cahn–Hilliard. Turbulent models implemented include RANS: Spalart-Allmaras and LES: Smagorinsky, Wale, Vreman; including wall models. Immersed boundary methods can be used, to avoid creating body fitted meshes. Acoustic propagation can be computed using Ffowcs-Williams and Hawkings models.
HORSES3D supports curvilinear, hexahedral, conforming meshes in GMSH, HDF5 and SpecMesh/HOHQMesh format. A hybrid CPU-based parallelisation strategy (shared and distributed memory) with OpenMP and MPI is followed.  

## External libraries
-----------
The following external routines/libraries can be used with **HORSES3D**, but are not necessary: METIS, MPI, HDF5, MKL, PETSc.  


## Compilers and third-party software
-----------
**HORSES3D** is an object-oriente Fortran 2008 solver, that can be compiled using gcc and intel compiler, in Unix-based operating systems. 
- We recommend using recent versions of such compilers (2019 or newer).
- Make is necessary (e.g., Gnu's version, which is included in most linux distributions).
- Supported meshes in GMSH, HDF5 and SpecMesh/HOHQMesh format.
- Post processing can be performed in tecplot or paraview.  


## Compiling & Running 
-----------
1. Go to the Solver folder to compile the code.  
        $ cd Solver  
2. Run configure file. This will configure the forders and test cases.  
        $ ./configure  
3. Install using the Make file (see manual in the /doc folder)  
        $ make clean  
        $ make all *Options*

	with the desired *Options* (default options are bold):  
         - MODE=DEBUG/**RELEASE**  
         - COMPILER=ifort/**gfortran**  
         - COMM=PARALLEL/**SEQUENTIAL**   
         - ENABLE THREADS=NO/**YES**  
         - WITH MKL=y/**n**  
    For example:  
         $ make all COMPILER=ifort COMM=PARALLEL  

4. Run the solver for the parameter file *file.control* (see manual in the /doc folder and examples in /test)  
        $ ./horses3d.ns *file.control*  
5. Test cases for various physics are provided in the folder /test  

## Additional libraries:  
-----------  
6. If you need to use HDF5, remember to set the paths  
        $ export HDF5_ROOT=path_to_hdf5  
        $ export HDF5_DIR=$HDF5_ROOT  
7. If you need to use PETSc, remember to set the path  
        $ export PETSC_DIR=path_to_petsc  
8. If you need to use METIS, remember to set the path  
        $ export METIS_HOME=path_to_metis  
