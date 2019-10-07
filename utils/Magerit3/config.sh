#!/bin/sh

module purge
module load petsc/3.11_Rintel2018b
module load HDF5/1.10.2-intel-2018b
export HDF5_ROOT=$HDF5_DIR
export METIS_HOME=/sw/software/METIS/5.1.0-foss-2018b

cd ~/wojtek/horses3d/Solver/
./configure
make all COMPILER=ifort COMM=PARALLEL
