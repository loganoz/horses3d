#################################################
#	Spectral Makefile 		#
#################################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J

LIBOBJS_NS =    ./build_ns/TimeIntegratorDefinitions.o \
		./build_ns/GMRESClass.o \
		./build_ns/JFNKClass.o \
		./build_ns/Implicit_JF.o \
		./build_ns/GenericLinSolverClass.o \
		./build_ns/PetscSolverClass.o \
		./build_ns/MKLPardisoSolverClass.o \
		./build_ns/IterativeSolverClass.o \
		./build_ns/MultigridTypes.o \
		./build_ns/TruncationErrorClass.o \
		./build_ns/ExplicitMethods.o \
		./build_ns/AnisFASMultigridClass.o \
		./build_ns/pAdaptationClass.o \
		./build_ns/FASMultigridClass.o \
		./build_ns/MultigridSolverClass.o \
		./build_ns/LinearSolverClass.o \
		./build_ns/AnalyticalJacobian.o \
		./build_ns/Implicit_NJ.o \
		./build_ns/IMEXMethods.o \
		./build_ns/TimeIntegrator.o 

LIBOBJS_CH =    ./build_ch/TimeIntegratorDefinitions.o \
		./build_ch/GMRESClass.o \
		./build_ch/JFNKClass.o \
		./build_ch/Implicit_JF.o \
		./build_ch/GenericLinSolverClass.o \
		./build_ch/PetscSolverClass.o \
		./build_ch/MKLPardisoSolverClass.o \
		./build_ch/IterativeSolverClass.o \
		./build_ch/MultigridTypes.o \
		./build_ch/TruncationErrorClass.o \
		./build_ch/ExplicitMethods.o \
		./build_ch/AnisFASMultigridClass.o \
		./build_ch/pAdaptationClass.o \
		./build_ch/FASMultigridClass.o \
		./build_ch/MultigridSolverClass.o \
		./build_ch/LinearSolverClass.o \
		./build_ch/AnalyticalJacobian.o \
		./build_ch/Implicit_NJ.o \
		./build_ch/IMEXMethods.o \
		./build_ch/TimeIntegrator.o 

LIB = timeintegrator

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:=all

all: lib$(LIB)_ns.a

./build_ns/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -I./include_ns $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ns -c $< -o $@ $(MOD_FLAG) ./include_ns

./build_ns/%.o: **/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -I./include_ns $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ns -c $< -o $@ $(MOD_FLAG) ./include_ns

./build_ch/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DCAHNHILLIARD -I./include_ch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ch -c $< -o $@ $(MOD_FLAG) ./include_ch

./build_ch/%.o: **/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DCAHNHILLIARD -I./include_ch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ch -c $< -o $@ $(MOD_FLAG) ./include_ch

lib$(LIB)_ns.a: header mkdirs $(LIBOBJS_NS) 
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_NS) 
	@$(RANLIB) $@

lib$(LIB)_ch.a: header mkdirs $(LIBOBJS_CH)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_CH) 
	@$(RANLIB) $@


install: all
	cp -p lib$(LIB)_ns.a $(INSTALL_DIR)/lib/
#	cp -p lib$(LIB)_ch.a $(INSTALL_DIR)/lib/
	cp -p ./include_ns/*.mod $(INSTALL_DIR)/include_ns/
#	cp -p ./include_ch/*.mod $(INSTALL_DIR)/include_ch/
	@echo
	@echo

header: FORCE
	@echo
	@echo "================================"
	@echo ">>   Building Time integrator    "
	@echo "================================"
	@echo
	@echo

mkdirs: FORCE
	@mkdir -p ./include_ns ./include_ch
	@mkdir -p ./build_ns ./build_ch

clean: FORCE
	$(RM) lib$(LIB)_ns.a  lib$(LIB)_ch.a
	$(RM) ./build_ns/*.o ./include_ns/*.mod
	$(RM) ./build_ch/*.o ./include_ch/*.mod

#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE:
