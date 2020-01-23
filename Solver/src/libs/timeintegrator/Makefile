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
		./build_ns/AnalyticalJacobian.o \
		./build_ns/GenericLinSolverClass.o \
		./build_ns/PetscSolverClass.o \
		./build_ns/GenericSmoother.o \
		./build_ns/MatrixFreeSmootherClass.o \
		./build_ns/MultigridSolverClass.o \
		./build_ns/MKLPardisoSolverClass.o \
		./build_ns/IterativeSolverClass.o \
		./build_ns/MatrixFreeGMRESClass.o \
		./build_ns/StaticCondensationSolverClass.o \
		./build_ns/LinearSolverClass.o \
		./build_ns/IMEXMethods.o \
		./build_ns/BDFTimeIntegrator.o \
		./build_ns/RosenbrockTimeIntegrator.o \
		./build_ns/MultigridTypes.o \
		./build_ns/TruncationErrorClass.o \
		./build_ns/ExplicitMethods.o \
		./build_ns/AnisFASMultigridClass.o \
		./build_ns/FASMultigridClass.o \
		./build_ns/MultiTauEstimationClass.o \
		./build_ns/pAdaptationClass.o \
		./build_ns/TimeIntegrator.o 

LIBOBJS_iNS =   ./build_ins/TimeIntegratorDefinitions.o \
		./build_ins/AnalyticalJacobian.o \
		./build_ins/GenericLinSolverClass.o \
		./build_ins/PetscSolverClass.o \
		./build_ins/GenericSmoother.o \
		./build_ins/MatrixFreeSmootherClass.o \
		./build_ins/MultigridSolverClass.o \
		./build_ins/MKLPardisoSolverClass.o \
		./build_ins/IterativeSolverClass.o \
		./build_ins/MatrixFreeGMRESClass.o \
		./build_ins/StaticCondensationSolverClass.o \
		./build_ins/LinearSolverClass.o \
		./build_ins/IMEXMethods.o \
		./build_ins/BDFTimeIntegrator.o \
		./build_ins/RosenbrockTimeIntegrator.o \
		./build_ins/MultigridTypes.o \
		./build_ins/TruncationErrorClass.o \
		./build_ins/ExplicitMethods.o \
		./build_ins/AnisFASMultigridClass.o \
		./build_ins/FASMultigridClass.o \
		./build_ins/MultiTauEstimationClass.o \
		./build_ins/pAdaptationClass.o \
		./build_ins/TimeIntegrator.o 


LIBOBJS_CH =    ./build_ch/TimeIntegratorDefinitions.o \
		./build_ch/AnalyticalJacobian.o \
		./build_ch/GenericLinSolverClass.o \
		./build_ch/PetscSolverClass.o \
		./build_ch/GenericSmoother.o \
		./build_ch/MatrixFreeSmootherClass.o \
		./build_ch/MultigridSolverClass.o \
		./build_ch/MKLPardisoSolverClass.o \
		./build_ch/IterativeSolverClass.o \
		./build_ch/MatrixFreeGMRESClass.o \
		./build_ch/StaticCondensationSolverClass.o \
		./build_ch/LinearSolverClass.o \
		./build_ch/IMEXMethods.o \
		./build_ch/BDFTimeIntegrator.o \
		./build_ch/RosenbrockTimeIntegrator.o \
		./build_ch/MultigridTypes.o \
		./build_ch/TruncationErrorClass.o \
		./build_ch/ExplicitMethods.o \
		./build_ch/AnisFASMultigridClass.o \
		./build_ch/FASMultigridClass.o \
		./build_ch/MultiTauEstimationClass.o \
		./build_ch/pAdaptationClass.o \
		./build_ch/TimeIntegrator.o 

LIBOBJS_MU =  ./build_mu/TimeIntegratorDefinitions.o \
		./build_mu/AnalyticalJacobian.o \
		./build_mu/GenericLinSolverClass.o \
		./build_mu/PetscSolverClass.o \
		./build_mu/GenericSmoother.o \
		./build_mu/MatrixFreeSmootherClass.o \
		./build_mu/MultigridSolverClass.o \
		./build_mu/MKLPardisoSolverClass.o \
		./build_mu/IterativeSolverClass.o \
		./build_mu/MatrixFreeGMRESClass.o \
		./build_mu/StaticCondensationSolverClass.o \
		./build_mu/LinearSolverClass.o \
		./build_mu/IMEXMethods.o \
		./build_mu/BDFTimeIntegrator.o \
		./build_mu/RosenbrockTimeIntegrator.o \
		./build_mu/MultigridTypes.o \
		./build_mu/TruncationErrorClass.o \
		./build_mu/ExplicitMethods.o \
		./build_mu/AnisFASMultigridClass.o \
		./build_mu/FASMultigridClass.o \
		./build_mu/MultiTauEstimationClass.o \
		./build_mu/pAdaptationClass.o \
		./build_mu/TimeIntegrator.o 

LIB = timeintegrator

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:=all

all: lib$(LIB)_ns.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_mu.a

./build_ns/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -I./include_ns $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ns -c $< -o $@ $(MOD_FLAG) ./include_ns

./build_ns/%.o: **/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -I./include_ns $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ns -c $< -o $@ $(MOD_FLAG) ./include_ns

./build_ins/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DINCNS -I./include_ins $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ins -c $< -o $@ $(MOD_FLAG) ./include_ins

./build_ins/%.o: **/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DINCNS -I./include_ins $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ins -c $< -o $@ $(MOD_FLAG) ./include_ins

./build_ch/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DCAHNHILLIARD -I./include_ch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ch -c $< -o $@ $(MOD_FLAG) ./include_ch

./build_ch/%.o: **/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DCAHNHILLIARD -I./include_ch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ch -c $< -o $@ $(MOD_FLAG) ./include_ch

./build_mu/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DMULTIPHASE -DCAHNHILLIARD -I./include_mu $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_mu -c $< -o $@ $(MOD_FLAG) ./include_mu

./build_mu/%.o: **/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DMULTIPHASE -DCAHNHILLIARD -I./include_mu $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_mu -c $<  -o $@ $(MOD_FLAG) ./include_mu

lib$(LIB)_ns.a: header mkdirs $(LIBOBJS_NS) 
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_NS) 
	@$(RANLIB) $@

lib$(LIB)_ins.a: header mkdirs $(LIBOBJS_iNS) 
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_iNS) 
	@$(RANLIB) $@

lib$(LIB)_ch.a: header mkdirs $(LIBOBJS_CH)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_CH) 
	@$(RANLIB) $@

lib$(LIB)_mu.a: header mkdirs $(LIBOBJS_MU)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_MU) 
	@$(RANLIB) $@

install: all
	cp -p lib$(LIB)_ns.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_ins.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_ch.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_mu.a $(INSTALL_DIR)/lib/
	cp -p ./include_ns/*.mod $(INSTALL_DIR)/include_ns/
	cp -p ./include_ins/*.mod $(INSTALL_DIR)/include_ins/
	cp -p ./include_ch/*.mod $(INSTALL_DIR)/include_ch/
	cp -p ./include_mu/*.mod $(INSTALL_DIR)/include_mu/
	@echo
	@echo

install-mu: lib$(LIB)_mu.a
	cp -p lib$(LIB)_mu.a $(INSTALL_DIR)/lib/
	cp -p ./include_mu/*.mod $(INSTALL_DIR)/include_mu/
	@echo
	@echo

install-ns: lib$(LIB)_ns.a
	cp -p lib$(LIB)_ns.a $(INSTALL_DIR)/lib/
	cp -p ./include_ns/*.mod $(INSTALL_DIR)/include_ns/
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
	@mkdir -p ./include_ns ./include_ins ./include_ch ./include_mu
	@mkdir -p ./build_ns ./build_ins ./build_ch ./build_mu

clean: FORCE
	$(RM) lib$(LIB)_ns.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_mu.a
	$(RM) ./build_ns/*.o ./include_ns/*.mod
	$(RM) ./build_ins/*.o ./include_ins/*.mod
	$(RM) ./build_ch/*.o ./include_ch/*.mod
	$(RM) ./build_mu/*.o ./include_mu/*.mod

#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE:
