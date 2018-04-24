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

LIBOBJS_NSCH =  ./build_nsch/TimeIntegratorDefinitions.o \
		./build_nsch/GMRESClass.o \
		./build_nsch/JFNKClass.o \
		./build_nsch/Implicit_JF.o \
		./build_nsch/GenericLinSolverClass.o \
		./build_nsch/PetscSolverClass.o \
		./build_nsch/MKLPardisoSolverClass.o \
		./build_nsch/IterativeSolverClass.o \
		./build_nsch/MultigridTypes.o \
		./build_nsch/TruncationErrorClass.o \
		./build_nsch/ExplicitMethods.o \
		./build_nsch/AnisFASMultigridClass.o \
		./build_nsch/pAdaptationClass.o \
		./build_nsch/FASMultigridClass.o \
		./build_nsch/MultigridSolverClass.o \
		./build_nsch/LinearSolverClass.o \
		./build_nsch/AnalyticalJacobian.o \
		./build_nsch/Implicit_NJ.o \
		./build_nsch/IMEXMethods.o \
		./build_nsch/TimeIntegrator.o 

LIB = timeintegrator

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:=all

all: lib$(LIB)_ns.a lib$(LIB)_ch.a lib$(LIB)_nsch.a

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

./build_nsch/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -DCAHNHILLIARD -I./include_nsch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nsch -c $< -o $@ $(MOD_FLAG) ./include_nsch

./build_nsch/%.o: **/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -DCAHNHILLIARD -I./include_nsch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nsch -c $<  -o $@ $(MOD_FLAG) ./include_nsch

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

lib$(LIB)_nsch.a: header mkdirs $(LIBOBJS_NSCH)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_NSCH) 
	@$(RANLIB) $@

install: all
	cp -p lib$(LIB)_ns.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_ch.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_nsch.a $(INSTALL_DIR)/lib/
	cp -p ./include_ns/*.mod $(INSTALL_DIR)/include_ns/
	cp -p ./include_ch/*.mod $(INSTALL_DIR)/include_ch/
	cp -p ./include_nsch/*.mod $(INSTALL_DIR)/include_nsch/
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
	@mkdir -p ./include_ns ./include_ch ./include_nsch
	@mkdir -p ./build_ns ./build_ch ./build_nsch

clean: FORCE
	$(RM) lib$(LIB)_ns.a  lib$(LIB)_ch.a lib$(LIB)_nsch.a
	$(RM) ./build_ns/*.o ./include_ns/*.mod
	$(RM) ./build_ch/*.o ./include_ch/*.mod
	$(RM) ./build_nsch/*.o ./include_nsch/*.mod

#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE:
