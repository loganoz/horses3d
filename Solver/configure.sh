#!/bin/bash

TEST_CASES="./Euler/BoxAroundCircle \
	    ./Euler/Diffuser \
	    ./Euler/PeriodicFlow \
	    ./Euler/UniformFlow \
	    ./Euler/UniformFlowPETSc \
	    ./Euler/JFNK \
            ./NavierStokes/Cylinder \
	    ./NavierStokes/FlatPlate \
	    ./NavierStokes/TaylorGreen \
	    ./NavierStokes/ManufacturedSolutions"

printf 'NSLITE3D_PATH = '$PWD'\nFTObject_PATH = '$PWD'/ftobjectlibrary' > ./test/make.inc

for subdir in $TEST_CASES; do
	mkdir -p ./test/$subdir/SETUP
	mkdir -p ./test/$subdir/RESULTS
	cp -v ./test/make.inc ./test/$subdir/make.inc
	cp -v ./test/Makefile.template ./test/$subdir/SETUP/Makefile
done


