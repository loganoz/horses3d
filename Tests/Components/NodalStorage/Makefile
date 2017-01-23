FC=gfortran
INC_DIR=./include
BUILD_DIR=./build
BIN_DIR=./bin
SRC_DIR=./src
EXECUTABLE = NodalStorage
MODE ?= DEBUG


OBJS_CASE_F90 = $(notdir $(wildcard $(SRC_DIR)/*.f90))
OBJS_CASE= $(patsubst %.f90,%,$(OBJS_CASE_F90))

OBJS_REMOVE = $(OBJS_CASE) NSLite3DMain ProblemFile TimeIntegrator 


# --------------- DO NOT EDIT ----------------------------------------------------

# --------------- CHOOSE OS -----------------------

include make.inc

# -------------- CHOOSE COMPILER ------------------
ifeq '$(FC)''gfortran'
   MOD_FLAG :=-J
   ifeq '$(MODE)''DEBUG'
      $(info MODE = DEBUG)
      #FFLAGS = -ffree-line-length-0 -O0 -g -Wall -fcheck=all -fbacktrace -fbounds-check -finit-real=snan -ffpe-trap=invalid,zero,overflow -fcheck-array-temporaries -fbackslash#-fdefault-integer-8 
      FFLAGS = -ffree-line-length-0 -cpp -O0 -g -fcheck=all -Wno-unused-variable -fbacktrace -fbounds-check -finit-real=snan  -ffpe-trap=invalid,zero,overflow -fcheck-array-temporaries -D_has_Quad
   else 
  	$(info MODE = RELEASE)
      FFLAGS = -ffree-line-length-0 -O3 -ftree-vectorize  -ftree-vectorizer-verbose=0 -fbackslash -fopenmp#-fdefault-integer-8 #-fopt-info-all
   endif
endif

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
	$(OBJS_EXEC)\

createFolders:
	mkdir -p $(BUILD_DIR)
	mkdir -p $(INC_DIR)
	mkdir -p $(BIN_DIR)


# ---------- Default compilation rule -----------
$(OBJS_CASE): %: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(INCLUDE) -J$(INC_DIR) -c $< -o $(BUILD_DIR)/$@.o
	


clean:
	@rm -f $(FTOBJLIB_PATH)/build/*
	@rm -f $(FTOBJLIB_PATH)/include/*
	@rm -f $(NSLITE3D_PATH)/build/*
	@rm -f $(NSLITE3D_PATH)/include/*
	@rm -f ./include/*.mod ./build/*.o
	@rm -f -r $(BIN_DIR)/*.dSYM
	
allclean: clean
	@rm -f $(BIN_DIR)/$(EXECUTABLE)

