#!/bin/bash

gfortran -o ascii_to_binary_vtk ascii_to_binary_vtk.f90

./ascii_to_binary_vtk ../PRECURSOR/17000/U_slice_boundary.vtk ../PRECURSOR/17000/U_slice_boundary_binary.vtkbin