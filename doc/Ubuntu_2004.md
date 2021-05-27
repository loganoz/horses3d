## Table of contents
* [General info](#general-info)
* [Installing dependencies](#dependencies)
* [Setup](#setup)
* [Known Issues](#issues)


## General info
This document is a quick guide to install HORSES3D in a vanilla enviroment with Ubuntu 20.04; as some extra steps may be required to add missing libraries and make Horse3D fully work. I'm not expert, and I don't control the code, so I may miss some extra things/libraries, or may be installing too much.

Note 0: only valid for Ubuntu 20.04. In previous versions versions some step were different, and some packages/libraries couldn't be compiled.

Note 1: In one setup I had some trouble with gfortran version, but I don't have it fully characterised. I just remmember having to switch back-and-forth from gfortran-7 to gfortran-9. Probably it was in my Ubuntu 18.04 set-up

Note 2: My makefile required some modifications. I added those steps to Makefile.in to remove manual modifications.

Note 3: Math libraries support:
If it's just for horse3D one can be chosen, depending on your processor, choose one. LAPACK is open-source?, MKL is from intel but may have better performance in intel processors.
Although HOPR can use MKL, I had some issues compiling it. I use lapack for HOPR and MKL in Horse3D

Note 4:
I lost some records on my .bash_history, so some extra dependencies may be required



## Installing dependencies
I may have installed to much, when trying to compile HORSE3D on my own. As this steps have worked in vanilla installations, I just keep all the steps.


### Basic tools and repo update

Installing the basic compiling tools, as used by OpenFOAM and others
```
sudo apt-get update; sudo apt-get upgrade
sudo apt-get install build-essential
sudo apt install bison flex m4
```

some extra build tools, required in some step (HOPR mainly):
```
sudo apt-get install cmake cmake-curses-gui
```


### basic dependencies for Horse3D:
--------
(don' use this)
~~sudo apt-get install gfortran gfortran-7 gfortran-7-multilib~~


### basic dependencies for HOPR:
Note that compilation requires python to be python2. Afterwards, switch back

```
sudo apt-get install gfortran-9 gfortran-9-multilib
sudo apt-get install liblapack-dev zlib1g-dev libc6 zlib* zlibc
sudo apt-get install python-is_python2
```

### HDF5 support for HOPR & HORSES3D
```
sudo apt-get install hdf5-tools hdf5-helpers
sudo apt-get install h5utils
```

¡WARNING, SOME NUMBERS MAY CHANGE!! take these as couple of examples
ubuntu 20.04:
```
sudo apt-get install libhdf5-103 libhdf5-cpp-103 libhdf5-dev
```
ubuntu 18.04
```
sudo apt-get install libhdf5-100 libhdf5-cpp-100 libhdf5-dev 
```



### MKL:
!WARNING not available in ubuntu 18.04
```
sudo apt-get install intel-mkl
```


### MPI:
```
sudo apt-get install openmpi-common openmpi-bin
sudo apt-get install libhdf5-mpi* libhdf5-openmpi-* libopenmpi-dev
```

### Forgot for what exactly
```
sudo apt-get install szip 
sudo apt-get install netcdf-bin curl zlib*
```


### PETSC:
```
sudo apt-get install petsc-dev
```


### METIS:
```
sudo apt-get install metis metis-edf libmetis-dev
```





## Setup


### HOPR (mesher):

Version 1 and 2 from the official website require some less work
```
make
```

Github version:
```
	mkdir build
	cd build/
	CC= FC= ccmake ../
		- c to configure (remove MPI, it failed)
		- c to configure
		- g to generate
	make
	../../bin/hopr parameter.ini
```

### Horse3D compile:

```
./configure
make clean
make all COMPILER=gfortran WITH_MKL=YES WITH_HDF5=YES COMM=PARALLEL WITH_PETSC=YES
```

add an alias to .bashrc, yo may be using it quite a few times (even compiling Problemfile.90)


## Known Issues
(TODO)
If using a h5 mesh fails with message "HD5F not link correctly"

* checking the path:
ldconfig -p | grep libhdf5.so
dpkg -L libhdf5-dev


* checking the path:
ldconfig -p | grep petsc

FOrtran segmentation fault, tipical mistakes:



# Extract a single solution

../../horses3d_stable/Solver/bin/horses2plt RESULTS/Cylinder_0000000000.hsol MESH/Cylinder_0000000000.hmesh --output-order=1 --output-basis=Homogeneous