- [General info](#general-info)
- [Installing dependencies](#installing-dependencies)
  - [Basic compiling tools and repo update](#basic-compiling-tools-and-repo-update)
  - [## basic dependencies for Horse3D:](#-basic-dependencies-for-horse3d)
  - [HDF5 support for HOPR & HORSES3D (synthetic wind field)](#hdf5-support-for-hopr--horses3d-synthetic-wind-field)
  - [MKL:](#mkl)
  - [MPI:](#mpi)
  - [Forgot for what exactly](#forgot-for-what-exactly)
  - [PETSC:](#petsc)
  - [METIS:](#metis)
- [Setup / Compilation](#setup--compilation)
  - [HOPR (mesher) - NOT NEEDED, we use the alternative gmsh:](#hopr-mesher---not-needed-we-use-the-alternative-gmsh)
  - [Horse3D compile:](#horse3d-compile)
- [Known Issues](#known-issues)
- [Extract a single solution](#extract-a-single-solution)

# General info
This document is a quick guide to install HORSES3D in a vanilla enviroment with Ubuntu 20.04; as some extra steps may be required to add missing libraries and make Horse3D fully work. I'm not expert, and I don't control the code, so I may miss some extra things/libraries, or may be installing too much (I copied a partial .bash_history to recover dependencies).

Note 0: This guide has been tested on Ubuntu 20.04 (native + WSL). In Ubuntu 18.04 some dependencies weren't available in default repositories.

Note 2: My makefile required some modifications to make hdf5 libraries work. I added those steps to Makefile.in to remove manual modifications.
        This modifications may only be valid for Ubuntu 20.04

Note 3: Math libraries support:
If it's just for horse3D one can be chosen, depending on your processor, choose one. LAPACK is open-source?, MKL is from intel but may have better performance in intel processors.
Although HOPR can use MKL, I had some issues compiling it. I use lapack for HOPR and MKL in Horse3D


# Installing dependencies
I may have installed to much, when trying to compile HORSE3D on my own. As this steps have worked in vanilla installations, I just keep all the steps.


## Basic compiling tools and repo update

Installing the basic compiling tools, as used by OpenFOAM and others
```bash
sudo apt-get update; sudo apt-get upgrade
sudo apt-get install build-essential
sudo apt install bison flex m4
```

some extra build tools, required in some step (HOPR mainly):
```bash
sudo apt-get install cmake cmake-curses-gui
```


## basic dependencies for Horse3D:
--------

```bash
sudo apt-get install gfortran-9 gfortran-9-multilib
sudo apt-get install liblapack-dev zlib1g-dev libc6 zlib* zlibc
```

## HDF5 support for HOPR & HORSES3D (synthetic wind field)
```bash
sudo apt-get install hdf5-tools hdf5-helpers
sudo apt-get install h5utils
```

¡WARNING, SOME NUMBERS MAY CHANGE!! take these as couple of examples
ubuntu 20.04:
```bash
sudo apt-get install libhdf5-103 libhdf5-cpp-103 libhdf5-dev
```


## MKL:
```bash
sudo apt-get install intel-mkl
```

Note: A warning will appear asking whether you want to use intel-MKL (propietary) as the default alternative to BLAS/LAPACK (open source). We have decided to say no when using AMD processor.

Note: See Note 3 in [General info](#general-info). Compiling HOPR with MKL may present issues. But can be manually overriden in its compilation.


## MPI:
To enable multithreading capabilities in HORSES3D
```bash
sudo apt-get install openmpi-common openmpi-bin
sudo apt-get install libhdf5-mpi* libhdf5-openmpi-* libopenmpi-dev
```

## Forgot for what exactly
```bash
sudo apt-get install netcdf-bin curl zlib*
```


## PETSC:
```bash
sudo apt-get install petsc-dev
```


## METIS:
```bash
sudo apt-get install metis metis-edf libmetis-dev
```




# Setup / Compilation

## HOPR (mesher) - NOT NEEDED, we use the alternative gmsh:

https://www.hopr-project.org

note: requires to install

```bash
sudo apt-get install python-is_python2
```

Version 1 and 2 from the official website require some less work
```bash
make
```

Github version:
```bash
	mkdir build
	cd build/
	CC= FC= ccmake ../
		- c to configure (remove MPI, it failed)
		- c to configure
		- g to generate
	make
	../../bin/hopr parameter.ini
```

## Horse3D compile:

Go to the folder Solver.

```bash
./configure
make clean
make all COMPILER=gfortran WITH_MKL=YES WITH_HDF5=YES COMM=PARALLEL WITH_PETSC=YES
```

If this gives an error:
-bash: ./configure: /bin/bash^M: bad interpreter: No such file or directory
run:
```bash
dos2unix configure
```
and try again. This error can also appear when compiling a makefile.


add an alias to .bashrc, yo may be using it quite a few times (even compiling Problemfile.90)


# Known Issues
(TODO)
If using a h5 mesh fails with message "HD5F not link correctly"

* checking the path:

```bash
ldconfig -p | grep libhdf5.so
dpkg -L libhdf5-dev

h5fc -show
```

* checking the path:
```bash
ldconfig -p | grep petsc
```

Fortran segmentation fault, tipical mistakes:



# Extract a single solution

../../horses3d_stable/Solver/bin/horses2plt RESULTS/Cylinder_0000000000.hsol MESH/Cylinder_0000000000.hmesh --output-order=1 --output-basis=Homogeneous