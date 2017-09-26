#!/bin/sh -f

FC=@CMAKE_Fortran_COMPILER@

cmd="$FC $* @CMAKE_Fortran_FLAGS@ @HORSES_F90FLAGS@ @CMAKE_SHARED_LIBRARY_Fortran_FLAGS@ @CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS@ -I@INCDIR@ -L@LIBDIR@ -lphysics -lftobject"
printf "%s " $cmd
printf "\n"
$FC $* @CMAKE_Fortran_FLAGS@ @HORSES_F90FLAGS@ @CMAKE_SHARED_LIBRARY_Fortran_FLAGS@ @CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS@ -I@INCDIR@ -L@LIBDIR@ -lphysics -lftobject

