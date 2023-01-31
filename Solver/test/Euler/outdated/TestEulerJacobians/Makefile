#####################################################################
## HORSES - a High-Order Spectral Element Solver
##
##      Test case makefile. This Makefile compiles a shared library,
##   libproblemfile.so with user specific procedures, such as initial
##   and boundary conditions or solution checks.
##
##  Usage:
##    -> cd $(CASE_DIR)/SETUP
##    -> make COMPILER=gfortran/ifort MODE=RELEASE/DEBUG COMM=SERIAL/PARALLEL
##    -> cd $(CASE_DIR)
##    -> ldd ./HORSES3D (to check whether the solver points
##        to the default library (found in
##        $HOME_DIR) or to the local library,
##        ./SETUP/libproblemfile.so)
##
#####################################################################
HOME_DIR=../../..
########################################
# 1.- User defined parameters
########################################
PLATFORM?=LINUX
COMPILER?=gfortran
MODE?=RELEASE
COMM?=SEQUENTIAL
ENABLE_THREADS?=YES
MAKE=make

########################################
# 2.- User defined paths
########################################
HORSESLIB=-L$(HOME_DIR)/lib -lphysicsns -lftobject -lio -lfoundation -lmpiutils
HORSESINC=$(HOME_DIR)/include
NSINC=$(HOME_DIR)/include_ns

########################################
# 3.- User defined compilation flags 
########################################
FFLAGS=-fPIC
GNU_RELEASE_FLAGS= -cpp -ffree-line-length-0 -O3 -ftree-vectorize  -ftree-vectorizer-verbose=0 \
                   -fbackslash
GNU_DEBUG_FLAGS=  -ffree-line-length-0 -cpp -O0 -g -fcheck=all -Wno-unused-variable \
                  -fbacktrace -fbounds-check -finit-real=snan -Wall \
                  -ffpe-trap=invalid,zero,overflow -fcheck-array-temporaries \
                  -fbackslash -fcheck=pointer
INTEL_RELEASE_FLAGS= -cpp -O3 -assume bscc -fp-model precise
INTEL_DEBUG_FLAGS= -cpp -O0 -g -warn all -fpscomp logicals -traceback -check all \
                   -check bounds -check uninit -ftrapuv -debug all -warn interfaces \
                   -ftrace=full -ftz -fp-model precise \
                   -fp-speculation=off -assume bscc
ARFLAGS=crv

TARGET=./libproblemfile

##########################################
## Select compiler 
##########################################
ifeq ($(COMPILER),gfortran)
  FC:=gfortran
  MOD_FLAG:=-J
  WHOLE_FILE:=--whole-file
else ifeq ($(COMPILER),gnu)
  FC:=gfortran
  MOD_FLAG:=-J
  WHOLE_FILE:=--whole-file
else ifeq ($(COMPILER),GNU)
  FC:=gfortran
  MOD_FLAG:=-J
  WHOLE_FILE:=--whole-file
endif

ifeq ($(COMPILER),ifort)
  FC:=ifort
  MOD_FLAG:=-module
  WHOLE_FILE:=
else ifeq ($(COMPILER),intel)
  FC:=ifort
  MOD_FLAG:=-module
  WHOLE_FILE:=
else ifeq ($(COMPILER),INTEL)
  FC:=ifort
  MOD_FLAG:=-module
  WHOLE_FILE:=
endif


all:
  $(FC) -I$(HORSESINC) -I$(NSINC) -o TestJacobians src/main.f90 $(HORSESLIB)
