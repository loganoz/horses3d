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
COMM?=SEQUENTIAL
ENABLE_THREADS?=YES
MAKE=make

########################################
# 2.- User defined paths
########################################
INCLUDE=-I$(HOME_DIR)/include
LIBS =
MACROS =
-include $(HOME_DIR)/make.inc
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
INTEL_RELEASE_FLAGS= -cpp -O3 -assume bscc -fp-model precise
INTEL_DEBUG_FLAGS= -cpp -O0 -g -warn all -fpscomp logicals -traceback -check all \
                   -check bounds -check uninit -ftrapuv -debug all -warn interfaces \
                   -ftrace=full -ftz -fp-model precise \
                   -fp-speculation=off -assume bscc
ARFLAGS=crv

TARGET=./libproblemfile

##########################################
## Select platform
##########################################
ifeq ($(PLATFORM),MACOSX)
   LIBSUFFIX=.dylib
   TARGET_NS=$(TARGET)_ns$(LIBSUFFIX)
   TARGET_CH=$(TARGET)_ch$(LIBSUFFIX)
   TARGET_NSCH=$(TARGET)_nsch$(LIBSUFFIX)
   DYNLIB_FLAG_NS=-fPIC -shared -install_name @rpath/$(TARGET_NS) -o $(TARGET_NS)
   DYNLIB_FLAG_CH=-fPIC -shared -install_name @rpath/$(TARGET_CH) -o $(TARGET_CH)
   DYNLIB_FLAG_NSCH=-fPIC -shared -install_name @rpath/$(TARGET_NSCH) -o $(TARGET_NSCH)
else ifeq ($(PLATFORM),LINUX)
   LIBSUFFIX=.so
   TARGET_NS=$(TARGET)_ns$(LIBSUFFIX)
   TARGET_CH=$(TARGET)_ch$(LIBSUFFIX)
   TARGET_NSCH=$(TARGET)_nsch$(LIBSUFFIX)
   DYNLIB_FLAG_NS=-fPIC -shared -o $(TARGET_NS)
   DYNLIB_FLAG_CH=-fPIC -shared -o $(TARGET_CH)
   DYNLIB_FLAG_NSCH=-fPIC -shared -o $(TARGET_NSCH)
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

ifeq ($(COMM),PARALLEL)
    ifeq ($(FC),gfortran)
        FC=mpif90
	MACROS+= -D_HAS_MPI_
    else ifeq ($(FC),ifort)
        FC=mpiifort
	MACROS+= -D_HAS_MPI_
    endif
endif

ifeq ($(ENABLE_THREADS),YES)
    FFLAGS+=-fopenmp
endif

##############################################################
##	Get objects 
##############################################################
HORSES_LIBS_NS=-lmonitors_ns -ljacobian_ns -lmesh_ns -lphysics_ns -lphysicsns -lio -lspectral -lmpiutils -lfoundation -lftobject
HORSES_LIBS_CH=-lmonitors_ch -ljacobian_ch -lmesh_ch -lphysics_ch -lphysicsch -lio -lspectral -lmpiutils -lfoundation -lftobject
HORSES_LIBS_NSCH=-lmonitors_ch -ljacobian_ch -lmesh_ch -lphysics_nsch -lphysicsch -lphysicsns -lio -lspectral -lmpiutils -lfoundation -lftobject

.DEFAULT_GOAL := all

all: $(TARGET_NS) $(TARGET_CH) $(TARGET_NSCH)

$(TARGET_NS): FORCE ./ProblemFile_NS.o
	$(FC) $(WHOLE_FILE) $(DYNLIB_FLAG_NS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_ns ./ProblemFile_NS.o -L$(HOME_DIR)/lib $(HORSES_LIBS_NS)

$(TARGET_CH): FORCE ./ProblemFile_CH.o
	$(FC) $(WHOLE_FILE) $(DYNLIB_FLAG_CH) $(FFLAGS) $(INCLUDE) $(INCLUDE)_ch ./ProblemFile_CH.o -L$(HOME_DIR)/lib $(HORSES_LIBS_CH)
	
$(TARGET_NSCH): FORCE ./ProblemFile_NSCH.o
	$(FC) $(WHOLE_FILE) $(DYNLIB_FLAG_NSCH) $(FFLAGS) $(INCLUDE) $(INCLUDE)_nsch ./ProblemFile_NSCH.o -L$(HOME_DIR)/lib $(HORSES_LIBS_NSCH)
	
# ---------- Executable macro -------------------
./ProblemFile_NS.o: FORCE ./ProblemFile.f90
	$(FC) -fPIC -shared -DNAVIERSTOKES $(FFLAGS) $(INCLUDE) $(INCLUDE)_ns -c ./ProblemFile.f90 -o $@

./ProblemFile_CH.o: FORCE ./ProblemFile.f90
	$(FC) -fPIC -shared -DCAHNHILLIARD $(FFLAGS) $(INCLUDE) $(INCLUDE)_ch -c ./ProblemFile.f90 -o $@

./ProblemFile_NSCH.o: FORCE ./ProblemFile.f90
	$(FC) -fPIC -shared -DNAVIERSTOKES -DCAHNHILLIARD $(FFLAGS) $(INCLUDE) $(INCLUDE)_nsch -c ./ProblemFile.f90 -o $@

clean: FORCE
	rm -f ./ProblemFile_NS.o ./ProblemFile_CH.o ./ProblemFile_NSCH.o ./*.mod
	rm -f ./*__genmod.f90
	
allclean: FORCE clean
	rm -f ./libproblemfile.so
	rm -f ./libproblemfile.dylib
	rm -f ./libproblemfile.dSYM
	rm -f ./libproblemfile_ns.so
	rm -f ./libproblemfile_ns.dylib
	rm -f ./libproblemfile_ns.dSYM
	rm -f ./libproblemfile_ch.so
	rm -f ./libproblemfile_ch.dylib
	rm -f ./libproblemfile_ch.dSYM
	rm -f ./libproblemfile_nsch.so
	rm -f ./libproblemfile_nsch.dylib
	rm -f ./libproblemfile_nsch.dSYM

FORCE:
