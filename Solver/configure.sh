#!/bin/bash

mkdir -p include
mkdir -p build
mkdir -p ftobjectlibrary
tar xf ftobjectlibrary.tar.gz -C ./ftobjectlibrary
mkdir -p ./test/Euler/BoxAroundCircle/RestartFiles
mkdir -p ./test/Euler/BoxAroundCircle/PlotFiles
mkdir -p ./test/Euler/diffuser/RestartFiles
mkdir -p ./test/Euler/diffuser/PlotFiles
mkdir -p ./test/Euler/PeriodicFlow/RestartFiles
mkdir -p ./test/Euler/PeriodicFlow/PlotFiles
mkdir -p ./test/Euler/UniformFlow/RestartFiles
mkdir -p ./test/Euler/UniformFlow/PlotFiles
mkdir -p ./test/Euler/UniformFlowPETSc/RestartFiles
mkdir -p ./test/Euler/UniformFlowPETSc/PlotFiles
mkdir -p ./test/Euler/JFNK/RestartFiles
mkdir -p ./test/Euler/JFNK/PlotFiles
mkdir -p ./test/NavierStokes/Cylinder/RestartFiles
mkdir -p ./test/NavierStokes/Cylinder/PlotFiles
mkdir -p ./test/NavierStokes/FlatPlate/RestartFiles
mkdir -p ./test/NavierStokes/FlatPlate/PlotFiles
mkdir -p ./test/NavierStokes/TaylorGreen/RestartFiles
mkdir -p ./test/NavierStokes/TaylorGreen/PlotFiles
mkdir -p ./test/NavierStokes/ManufacturedSolutions/RestartFiles
mkdir -p ./test/NavierStokes/ManufacturedSolutions/PlotFiles
printf 'NSLITE3D_PATH = '$PWD'\nFTObject_PATH = '$PWD'/ftobjectlibrary' > ./test/make.inc
cp -v ./test/make.inc ./test/Components/FacePatches/make.inc
cp -v ./test/make.inc ./test/Components/Gradients/make.inc
cp -v ./test/make.inc ./test/Components/HexMappings/make.inc
cp -v ./test/make.inc ./test/Components/HexMesh/make.inc
cp -v ./test/make.inc ./test/Components/MappedGeometry/make.inc
cp -v ./test/make.inc ./test/Components/NodalStorage/make.inc
cp -v ./test/make.inc ./test/Euler/BoxAroundCircle/make.inc
cp -v ./test/make.inc ./test/Euler/diffuser/make.inc
cp -v ./test/make.inc ./test/Euler/PeriodicFlow/make.inc
cp -v ./test/make.inc ./test/Euler/UniformFlow/make.inc
cp -v ./test/make.inc ./test/Euler/UniformFlowPETSc/make.inc
cp -v ./test/make.inc ./test/Euler/JFNK/make.inc
cp -v ./test/make.inc ./test/NavierStokes/Cylinder/make.inc
cp -v ./test/make.inc ./test/NavierStokes/FlatPlate/make.inc
cp -v ./test/make.inc ./test/NavierStokes/TaylorGreen/make.inc
cp -v ./test/make.inc ./test/NavierStokes/ManufacturedSolutions/make.inc

cp -v ./test/Makefile.template ./test/Euler/BoxAroundCircle/Makefile
cp -v ./test/Makefile.template ./test/Euler/diffuser/Makefile
cp -v ./test/Makefile.template ./test/Euler/PeriodicFlow/Makefile
cp -v ./test/Makefile.template ./test/Euler/UniformFlow/Makefile
cp -v ./test/Makefile.template ./test/Euler/UniformFlowPETSc/Makefile
cp -v ./test/Makefile.template ./test/Euler/JFNK/Makefile
cp -v ./test/Makefile.template ./test/NavierStokes/Cylinder/Makefile
cp -v ./test/Makefile.template ./test/NavierStokes/FlatPlate/Makefile
cp -v ./test/Makefile.template ./test/NavierStokes/TaylorGreen/Makefile
cp -v ./test/Makefile.template ./test/NavierStokes/ManufacturedSolutions/Makefile
