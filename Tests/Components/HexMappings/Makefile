COMPILER ?= gfortran
MODE ?= RELEASE

INC_DIR=./include
BUILD_DIR=./build
BIN_DIR=./bin
SRC_DIR=./src
EXECUTABLE = HexMappings


OBJS_CASE_F90 = $(notdir $(wildcard $(SRC_DIR)/*.f90))
OBJS_CASE= $(patsubst %.f90,%,$(OBJS_CASE_F90))

OBJS_REMOVE = $(OBJS_CASE) NSLite3DMain ProblemFile TimeIntegrator 


# --------------- DO NOT EDIT ----------------------------------------------------

# --------------- CHOOSE OS -----------------------

include make.inc

# -------------- CHOOSE COMPILER ------------------
ifeq ($(COMPILER),gfortran)
   MOD_FLAG :=-J
   FC=gfortran
   ifeq ($(MODE),DEBUG)
      $(info Compiling with gfortran in DEBUG mode.)
      FFLAGS = -ffree-line-length-0 -O0 -g -fcheck=all -Wno-unused-variable -fbacktrace -fbounds-check -finit-real=snan  -ffpe-trap=invalid,zero,overflow -fcheck-array-temporaries
   else 
     $(info Compiling with gfortran in RELEASE mode.)
     FFLAGS = -ffree-line-length-0 -O3 -ftree-vectorize  -ftree-vectorizer-verbose=0 -fbackslash -fopenmp#-fdefault-integer-8 #-fopt-info-all
   endif

else ifeq ($(COMPILER),ifort)
   MOD_FLAG := -module 
   FC=ifort
   ifeq ($(MODE),DEBUG)
     $(info Compiling with intel fortran in DEBUG mode.)
     FFLAGS = -O0 -g -warn all -fpscomp logicals -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interfaces -warn interfaces -ftrace=full -assume ieee_fpe_flags -ftz -fp-model precise -fp-speculation=off -assume bscc 
   else
     $(info Compiling with intel fortran in RELEASE mode.)
     FFLAGS = -O3 -assume bscc
   endif
endif
FFLAGS += -cpp -D_has_Quad


# ------------- LIBRARIES ------------------------
INCLUDE=-I$(INC_DIR)


# ------------- INVOKE NSLITE3D Objects --------------
include $(NSLITE3D_PATH)/Makefile.nslite3d

OBJS_CASE_EXEC = $(foreach obj,$(OBJS_CASE),./build/$(obj).o) 


# ---------- Remove the replaced files ---------------------

OBJS_NSLITE3D_Updated = $(filter-out $(OBJS_REMOVE),$(OBJS_NSLITE3D))
OBJS_NSLITE3D_EXEC_Updated = $(foreach obj,$(OBJS_NSLITE3D_Updated),$(NSLITE3D_PATH)/build/$(obj).o)

# -----------------------------------------------------------

OBJS = $(OBJS_FTOBJLIB) $(OBJS_NSLITE3D_Updated) $(OBJS_CASE)
OBJS_EXEC = $(OBJS_FTOBJLIB_EXEC) $(OBJS_NSLITE3D_EXEC_Updated) $(OBJS_CASE_EXEC)
	

all: createFolders create_$(EXECUTABLE)
	
# ---------- Executable macro -------------------
create_$(EXECUTABLE): $(OBJS) 
	$(FC) $(FFLAGS) -o $(BIN_DIR)/$(EXECUTABLE) $(INCLUDE) \
	$(OBJS_EXEC) $(LIBS)\

createFolders:
	mkdir -p $(BUILD_DIR)
	mkdir -p $(INC_DIR)
	mkdir -p $(BIN_DIR)



# ---------- Default compilation rule -----------
$(OBJS_CASE): %: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(INCLUDE) $(MOD_FLAG)$(INC_DIR) -c $< -o $(BUILD_DIR)/$@.o
	


clean:
	@rm -f $(FTObject_PATH)/build/*
	@rm -f $(FTObject_PATH)/include/*
	@rm -f $(NSLITE3D_PATH)/build/*
	@rm -f $(NSLITE3D_PATH)/include/*
	@rm -f ./include/*.mod ./build/*.o
	@rm -f -r $(BIN_DIR)/*.dSYM
	
allclean: clean
	@rm -f $(BIN_DIR)/$(EXECUTABLE)

