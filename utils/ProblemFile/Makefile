#####################################################################
## HORSES - a High-Order Spectral Element Solver
##
##      Test case makefile. This Makefile compiles a shared library,
##   libproblemfile.so with user specific procedures, such as initial
##   and boundary conditions or solution checks.
##
##	Usage:
##		-> cd $(CASE_DIR)/SETUP
##		-> make COMPILER=gfortran/ifort MODE=RELEASE/DEBUG COMM=SERIAL/PARALLEL
##		-> cd $(CASE_DIR)
##		-> ldd ./HORSES3D (to check whether the solver points
##				to the default library (found in
##				$HOME_DIR) or to the local library,
##				./SETUP/libproblemfile.so)
##
#####################################################################
HOME_DIR=@HORSES_HOME_DIR
########################################
# 1.- User defined parameters
########################################
PLATFORM?=LINUX
COMPILER?=gfortran
MODE?=RELEASE
COMM?=PARALLEL
MAKE=make

########################################
# 2.- User defined paths
########################################
INCLUDE=-I$(HOME_DIR)/include
LIBS =
MACROS =
-include make.inc

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
INTEL_RELEASE_FLAGS= -cpp -O3 -assume bscc
INTEL_DEBUG_FLAGS= -cpp -O0 -g -warn all -fpscomp logicals -traceback -check all \
                   -check bounds -check uninit -ftrapuv -debug all -warn interfaces \
                   -ftrace=full -assume ieee_fpe_flags -ftz -fp-model precise \
                   -fp-speculation=off -assume bscc
ARFLAGS=crv

##########################################
## Select platform
##########################################
ifeq ($(PLATFORM),LINUX)
     TARGET=./libproblemfile.so
     DYNLIB_FLAG=-fPIC -shared -o $(TARGET)
else ifeq ($(PLATFORM),MACOSX)
     TARGET=./libproblemfile.dylib
     DYNLIB_FLAG=-fPIC -shared -install_name @rpath/libproblemfile.dylib -o $(TARGET)
endif

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

##########################################
## Select mode
##########################################
ifeq ($(MODE),RELEASE)
     ifeq ($(FC),gfortran)
         FFLAGS+=$(GNU_RELEASE_FLAGS)
     else ifeq ($(FC),ifort)
         FFLAGS+=$(INTEL_RELEASE_FLAGS)
     endif
else ifeq ($(MODE),DEBUG)
     ifeq ($(FC),gfortran)
         FFLAGS+=$(GNU_DEBUG_FLAGS)
     else ifeq ($(FC),ifort)
         FFLAGS+=$(INTEL_DEBUG_FLAGS)
     endif
endif

ifeq ($(MODE),RELEASE)
   ifeq ($(COMM),PARALLEL)
      FFLAGS+=-fopenmp
   endif
endif


##############################################################
##	Get objects 
##############################################################
HORSES_LIBS=-lmonitors -ljacobian -lplot -lmesh -lphysics -lio -lspectral -lmpiutils -lfoundation -lftobject

.DEFAULT_GOAL := all

all: $(TARGET)

$(TARGET): ./ProblemFile.o
	$(FC) $(WHOLE_FILE) $(DYNLIB_FLAG) $(FFLAGS) $(INCLUDE) ./ProblemFile.o -L$(HOME_DIR)/lib $(HORSES_LIBS)
	
# ---------- Executable macro -------------------
./ProblemFile.o: ./ProblemFile.f90
	$(FC) -fPIC -shared $(FFLAGS) $(INCLUDE) -c ./ProblemFile.f90 -o $@

clean: FORCE
	rm -f ./ProblemFile.o ./*.mod
	rm -f ./*__genmod.f90
	
allclean: FORCE clean
	rm -f ./libproblemfile.so

FORCE:
