#####################################################################
## HORSES - a High-Order Spectral Element Solver
##
##	Main code Makefile
##
##
#####################################################################

########################################
# 1.- User defined parameters
########################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIN?=ranlib
MOD_FLAG?=-J
PLATFORM?=LINUX
COMPILER?=gfortran
MODE?=RELEASE
COMM?=PARALLEL
MAKE=make
AR=ar
RANLIB=ranlib
PROG=./bin/HORSES3D

OBJS_HORSES=	DGInviscidDiscretization \
		DGViscousDiscretization \
		DGWeakIntegrals \
		DGSEMClass \
		SpatialDiscretization \
		TimeIntegrator \
		pAdaptationClass \
		Implicit_JF \
		JFNKClass \
		GMRESClass \
		Implicit_NJ \
		PetscSolverClass \
		GenericLinSolverClass \
		MKLPardisoSolverClass \
		LinearSolverClass \
		IterativeSolverClass \
		MultigridSolverClass \
		ExplicitMethods \
		FASMultigridClass \
		HORSES3DMain 

OBJS_HORSES_BUILD=$(foreach obj,$(OBJS_HORSES),./build/$(obj).o)


##########################################
##	Include dependencies 		##
##########################################
-include ./make.deps


##########################################################
##		COMPILATION RULES 			##
##########################################################

.DEFAULT_GOAL := $(PROG)

./build/%.o: ./%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(INCLUDE) $(MACROS) -c $< -o $@ $(MOD_FLAG) ./include

$(PROG): header mkdirs $(OBJS_HORSES_BUILD)
	@echo
	@echo
	@echo "---------------------------------"
	@echo "| Linking " $@ "|"
	@echo "---------------------------------"
	$(FC) $(FFLAGS) -o $(PROG) $(OBJS_HORSES_BUILD) $(HORSES_LIBS_EXEC)


install: $(PROG)
	cp -p $(PROG) $(INSTALL_DIR)/bin
	cp -p ./include/*.mod $(INSTALL_DIR)/include
	@echo
	@echo

##########################################
##	Extra auxiliar rules 		##
##########################################
clean: FORCE
	rm -f ./include/*.mod ./lib/*.a ./build/*.o
	rm -r $(PROG)

mkdirs: FORCE
	@mkdir -p ./build
	@mkdir -p ./include
	@mkdir -p ./bin
	@mkdir -p ./lib

header: FORCE
	@echo
	@echo "================================"
	@echo ">>   Building HORSES3D Solver  "
	@echo "================================"
	@echo
	@echo

## Dummy target to make sure a rule is executed
FORCE:
