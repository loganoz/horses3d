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
   TARGET_iNS=$(TARGET)_ins$(LIBSUFFIX)
   TARGET_CH=$(TARGET)_ch$(LIBSUFFIX)
   TARGET_NSCH=$(TARGET)_nsch$(LIBSUFFIX)
   TARGET_iNSCH=$(TARGET)_insch$(LIBSUFFIX)
   DYNLIB_FLAG_NS=-fPIC -shared -install_name @rpath/$(TARGET_NS) -o $(TARGET_NS)
   DYNLIB_FLAG_iNS=-fPIC -shared -install_name @rpath/$(TARGET_iNS) -o $(TARGET_iNS)
   DYNLIB_FLAG_CH=-fPIC -shared -install_name @rpath/$(TARGET_CH) -o $(TARGET_CH)
   DYNLIB_FLAG_NSCH=-fPIC -shared -install_name @rpath/$(TARGET_NSCH) -o $(TARGET_NSCH)
   DYNLIB_FLAG_iNSCH=-fPIC -shared -install_name @rpath/$(TARGET_iNSCH) -o $(TARGET_iNSCH)
else ifeq ($(PLATFORM),LINUX)
   LIBSUFFIX=.so
   TARGET_NS=$(TARGET)_ns$(LIBSUFFIX)
   TARGET_iNS=$(TARGET)_ins$(LIBSUFFIX)
   TARGET_CH=$(TARGET)_ch$(LIBSUFFIX)
   TARGET_NSCH=$(TARGET)_nsch$(LIBSUFFIX)
   TARGET_iNSCH=$(TARGET)_insch$(LIBSUFFIX)
   DYNLIB_FLAG_NS=-fPIC -shared -o $(TARGET_NS)
   DYNLIB_FLAG_iNS=-fPIC -shared -o $(TARGET_iNS)
   DYNLIB_FLAG_CH=-fPIC -shared -o $(TARGET_CH)
   DYNLIB_FLAG_NSCH=-fPIC -shared -o $(TARGET_NSCH)
   DYNLIB_FLAG_iNSCH=-fPIC -shared -o $(TARGET_iNSCH)
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
HORSES_LIBS_NS=-ltimeintegrator_ns -ljacobian_ns -ldiscretization_ns -lmonitors_ns -lparticles_ns -lmesh_ns -lphysics_ns -lphysicsns -lio -lspectral -lmpiutils -lfoundation -lftobject
HORSES_LIBS_iNS=-ltimeintegrator_ins -ljacobian_ins -ldiscretization_ins -lmonitors_ins -lparticles_ins -lmesh_ins -lphysics_ins -lphysicsns -lio -lspectral -lmpiutils -lfoundation -lftobject
HORSES_LIBS_CH=-ltimeintegrator_ch -ljacobian_ch -ldiscretization_ch -lmonitors_ch -lparticles_ch -lmesh_ch -lphysics_ch -lphysicsch -lio -lspectral -lmpiutils -lfoundation -lftobject
HORSES_LIBS_NSCH=-ltimeintegrator_nsch -ljacobian_nsch -ldiscretization_nsch -lmonitors_nsch -lparticles_nsch -lmesh_nsch -lphysics_nsch -lphysicsch -lphysicsns -lio -lspectral -lmpiutils -lfoundation -lftobject
HORSES_LIBS_iNSCH=-ltimeintegrator_insch -ljacobian_insch -ldiscretization_insch -lmonitors_insch -lparticles_insch -lmesh_insch -lphysics_insch -lphysicsch -lphysicsins -lio -lspectral -lmpiutils -lfoundation -lftobject

.DEFAULT_GOAL := all

all: $(TARGET_NS) $(TARGET_iNS) $(TARGET_CH) $(TARGET_NSCH) $(TARGET_iNSCH)

$(TARGET_NS): FORCE ./ProblemFile_NS.o
	$(FC) $(WHOLE_FILE) $(DYNLIB_FLAG_NS) $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_ns ./ProblemFile_NS.o -L$(HOME_DIR)/lib $(HORSES_LIBS_NS)

$(TARGET_iNS): FORCE ./ProblemFile_iNS.o
	$(FC) $(WHOLE_FILE) $(DYNLIB_FLAG_iNS) $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_ins ./ProblemFile_iNS.o -L$(HOME_DIR)/lib $(HORSES_LIBS_iNS)

$(TARGET_CH): FORCE ./ProblemFile_CH.o
	$(FC) $(WHOLE_FILE) $(DYNLIB_FLAG_CH) $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_ch ./ProblemFile_CH.o -L$(HOME_DIR)/lib $(HORSES_LIBS_CH)
	
$(TARGET_NSCH): FORCE ./ProblemFile_NSCH.o
	$(FC) $(WHOLE_FILE) $(DYNLIB_FLAG_NSCH) $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_nsch ./ProblemFile_NSCH.o -L$(HOME_DIR)/lib $(HORSES_LIBS_NSCH)
	
$(TARGET_iNSCH): FORCE ./ProblemFile_iNSCH.o
	$(FC) $(WHOLE_FILE) $(DYNLIB_FLAG_iNSCH) $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_insch ./ProblemFile_iNSCH.o -L$(HOME_DIR)/lib $(HORSES_LIBS_iNSCH)
	
# ---------- Executable macro -------------------
./ProblemFile_NS.o: FORCE ./ProblemFile.f90
	$(FC) -fPIC -shared -DNAVIERSTOKES $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_ns -c ./ProblemFile.f90 -o $@

./ProblemFile_iNS.o: FORCE ./ProblemFile.f90
	$(FC) -fPIC -shared -DINCNS $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_ins -c ./ProblemFile.f90 -o $@

./ProblemFile_CH.o: FORCE ./ProblemFile.f90
	$(FC) -fPIC -shared -DCAHNHILLIARD $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_ch -c ./ProblemFile.f90 -o $@

./ProblemFile_NSCH.o: FORCE ./ProblemFile.f90
	$(FC) -fPIC -shared -DNAVIERSTOKES -DCAHNHILLIARD $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_nsch -c ./ProblemFile.f90 -o $@

./ProblemFile_iNSCH.o: FORCE ./ProblemFile.f90
	$(FC) -fPIC -shared -DINCNS -DCAHNHILLIARD $(MACROS) $(FFLAGS) $(INCLUDE) $(INCLUDE)_insch -c ./ProblemFile.f90 -o $@

clean: FORCE
	rm -f ./ProblemFile_NS.o ./ProblemFile_iNS.o ./ProblemFile_CH.o ./ProblemFile_NSCH.o ./ProblemFile_iNSCH.o ./*.mod
	rm -f ./*__genmod.f90
	
allclean: FORCE clean
	rm -f ./libproblemfile.so
	rm -f ./libproblemfile.dylib
	rm -f ./libproblemfile.dSYM
	rm -f ./libproblemfile_ns.so
	rm -f ./libproblemfile_ns.dylib
	rm -f ./libproblemfile_ns.dSYM
	rm -f ./libproblemfile_ins.so
	rm -f ./libproblemfile_ins.dylib
	rm -f ./libproblemfile_ins.dSYM
	rm -f ./libproblemfile_ch.so
	rm -f ./libproblemfile_ch.dylib
	rm -f ./libproblemfile_ch.dSYM
	rm -f ./libproblemfile_nsch.so
	rm -f ./libproblemfile_nsch.dylib
	rm -f ./libproblemfile_nsch.dSYM
	rm -f ./libproblemfile_insch.so
	rm -f ./libproblemfile_insch.dylib
	rm -f ./libproblemfile_insch.dSYM

FORCE:
