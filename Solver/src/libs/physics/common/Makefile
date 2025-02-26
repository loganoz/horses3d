#################################################
#	Foundation Makefile 		#
#################################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J

LIBOBJS_NS = 	./build_ns/FluidData.o \
		./build_ns/PhysicsStorage.o \
		./build_ns/VariableConversion.o \
		./build_ns/Physics.o \
		./build_ns/GenericBoundaryConditionClass.o \
		./build_ns/InflowBC.o \
		./build_ns/OutflowBC.o \
		./build_ns/NoSlipWallBC.o \
		./build_ns/FreeSlipWallBC.o \
		./build_ns/PeriodicBC.o \
		./build_ns/UserDefinedBC.o \
		./build_ns/LESModels.o \
		./build_ns/WallFunctionDefinitions.o \
		./build_ns/WallFunctionBC.o \
		./build_ns/ZYDIYFreeSlipWallBC.o \
		./build_ns/ZYDIYDirichletBC.o \
		./build_ns/BoundaryConditions.o 

LIBOBJS_NSSA = 	./build_nssa/FluidData.o \
		./build_nssa/PhysicsStorage.o \
		./build_nssa/VariableConversion.o \
		./build_nssa/Physics.o \
		./build_nssa/GenericBoundaryConditionClass.o \
		./build_nssa/InflowBC.o \
		./build_nssa/OutflowBC.o \
		./build_nssa/NoSlipWallBC.o \
		./build_nssa/FreeSlipWallBC.o \
		./build_nssa/PeriodicBC.o \
		./build_nssa/UserDefinedBC.o \
		./build_nssa/LESModels.o \
		./build_nssa/WallFunctionDefinitions.o \
		./build_nssa/WallFunctionBC.o \
		./build_nssa/ZYDIYFreeSlipWallBC.o \
		./build_nssa/ZYDIYDirichletBC.o \
		./build_nssa/BoundaryConditions.o 
		
LIBOBJS_iNS = 	./build_ins/FluidData.o \
		./build_ins/PhysicsStorage.o \
		./build_ins/VariableConversion.o \
		./build_ins/Physics.o \
		./build_ins/GenericBoundaryConditionClass.o \
		./build_ins/InflowBC.o \
		./build_ins/OutflowBC.o \
		./build_ins/NoSlipWallBC.o \
		./build_ins/FreeSlipWallBC.o \
		./build_ins/PeriodicBC.o \
		./build_ins/UserDefinedBC.o \
		./build_ins/LESModels.o \
		./build_ins/ZYDIYFreeSlipWallBC.o \
		./build_ins/ZYDIYDirichletBC.o \
		./build_ins/BoundaryConditions.o 

LIBOBJS_CH = 	./build_ch/FluidData.o \
		./build_ch/PhysicsStorage.o \
		./build_ch/VariableConversion.o \
		./build_ch/Physics.o \
		./build_ch/GenericBoundaryConditionClass.o \
		./build_ch/InflowBC.o \
		./build_ch/OutflowBC.o \
		./build_ch/NoSlipWallBC.o \
		./build_ch/FreeSlipWallBC.o \
		./build_ch/PeriodicBC.o \
		./build_ch/UserDefinedBC.o \
		./build_ch/ZYDIYFreeSlipWallBC.o \
		./build_ch/ZYDIYDirichletBC.o \
		./build_ch/BoundaryConditions.o 

LIBOBJS_MU = 	./build_mu/FluidData.o \
		./build_mu/PhysicsStorage.o \
		./build_mu/VariableConversion.o \
		./build_mu/Physics.o \
		./build_mu/GenericBoundaryConditionClass.o \
		./build_mu/InflowBC.o \
		./build_mu/OutflowBC.o \
		./build_mu/NoSlipWallBC.o \
		./build_mu/FreeSlipWallBC.o \
		./build_mu/PeriodicBC.o \
		./build_mu/UserDefinedBC.o \
		./build_mu/ZYDIYFreeSlipWallBC.o \
		./build_mu/ZYDIYDirichletBC.o \
		./build_mu/BoundaryConditions.o 

LIBOBJS_SLR = 	./build_slr/FluidData.o \
		./build_slr/PhysicsStorage.o \
		./build_slr/VariableConversion.o \
		./build_slr/Physics.o \
		./build_slr/GenericBoundaryConditionClass.o \
		./build_slr/InflowBC.o \
		./build_slr/OutflowBC.o \
		./build_slr/NoSlipWallBC.o \
		./build_slr/FreeSlipWallBC.o \
		./build_slr/PeriodicBC.o \
		./build_slr/UserDefinedBC.o \
		./build_slr/ZYDIYFreeSlipWallBC.o \
		./build_slr/ZYDIYDirichletBC.o \
		./build_slr/BoundaryConditions.o 

LIBOBJS_SLR_INS_V04 = 	./build_slr_ins_v04/FluidData.o \
		./build_slr_ins_v04/PhysicsStorage.o \
		./build_slr_ins_v04/VariableConversion.o \
		./build_slr_ins_v04/Physics.o \
		./build_slr_ins_v04/GenericBoundaryConditionClass.o \
		./build_slr_ins_v04/InflowBC.o \
		./build_slr_ins_v04/OutflowBC.o \
		./build_slr_ins_v04/NoSlipWallBC.o \
		./build_slr_ins_v04/FreeSlipWallBC.o \
		./build_slr_ins_v04/PeriodicBC.o \
		./build_slr_ins_v04/UserDefinedBC.o \
		./build_slr_ins_v04/ZYDIYFreeSlipWallBC.o \
		./build_slr_ins_v04/ZYDIYDirichletBC.o \
		./build_slr_ins_v04/BoundaryConditions.o 

LIB = physics

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL := lib$(LIB).a

all: lib$(LIB)_ns.a lib$(LIB)_nssa.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_mu.a   lib$(LIB)_slr_ins_v04.a

./build_ns/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -I./include_ns $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ns -c $< -o $@ $(MOD_FLAG) ./include_ns

./build_nssa/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -DSPALARTALMARAS -I./include_nssa $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nssa -c $< -o $@ $(MOD_FLAG) ./include_nssa

./build_ins/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DINCNS -I./include_ins $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ins -c $< -o $@ $(MOD_FLAG) ./include_ins

./build_ch/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DCAHNHILLIARD -I./include_ch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ch -c $< -o $@ $(MOD_FLAG) ./include_ch

./build_mu/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DMULTIPHASE -DCAHNHILLIARD -I./include_mu $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_mu -c $< -o $@ $(MOD_FLAG) ./include_mu

./build_slr/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DSCALAR -I./include_slr $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_slr -c $< -o $@ $(MOD_FLAG) ./include_slr


./build_slr_ins_v04/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DSCALAR_INS_V04 -I./include_slr_ins_v04 $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_slr_ins_v04 -c $< -o $@ $(MOD_FLAG) ./include_slr_ins_v04


lib$(LIB)_ns.a: header mkdirs $(LIBOBJS_NS) 
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_NS) 
	@$(RANLIB) $@

lib$(LIB)_nssa.a: header mkdirs $(LIBOBJS_NSSA) 
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_NSSA) 
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

lib$(LIB)_slr.a: header mkdirs $(LIBOBJS_SLR)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_SLR) 
	@$(RANLIB) $@

lib$(LIB)_slr_ins_v04.a: header mkdirs $(LIBOBJS_SLR_INS_V04)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_SLR_INS_V04) 
	@$(RANLIB) $@

install: all
	cp -p lib$(LIB)_ns.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_nssa.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_ins.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_ch.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_mu.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_slr_ins_v04.a $(INSTALL_DIR)/lib/
	cp -p ./include_ns/*.mod $(INSTALL_DIR)/include_ns/
	cp -p ./include_nssa/*.mod $(INSTALL_DIR)/include_nssa/
	cp -p ./include_ins/*.mod $(INSTALL_DIR)/include_ins/
	cp -p ./include_ch/*.mod $(INSTALL_DIR)/include_ch/
	cp -p ./include_mu/*.mod $(INSTALL_DIR)/include_mu/
	cp -p ./include_slr_ins_v04/*.mod $(INSTALL_DIR)/include_slr_ins_v04/
	@echo
	@echo

install-mu: lib$(LIB)_mu.a
	cp -p lib$(LIB)_mu.a $(INSTALL_DIR)/lib/
	cp -p ./include_mu/*.mod $(INSTALL_DIR)/include_mu/
	@echo
	@echo

install-slr: lib$(LIB)_slr.a
	cp -p lib$(LIB)_slr.a $(INSTALL_DIR)/lib/
	cp -p ./include_slr/*.mod $(INSTALL_DIR)/include_slr/
	@echo
	@echo

install-slr_ins_v04: lib$(LIB)_slr_ins_v04.a
	cp -p lib$(LIB)_slr_ins_v04.a $(INSTALL_DIR)/lib/
	cp -p ./include_slr_ins_v04/*.mod $(INSTALL_DIR)/include_slr_ins_v04/
	@echo
	@echo

install-ns: lib$(LIB)_ns.a
	cp -p lib$(LIB)_ns.a $(INSTALL_DIR)/lib/
	cp -p ./include_ns/*.mod $(INSTALL_DIR)/include_ns/
	@echo
	@echo

install-nssa: lib$(LIB)_nssa.a
	cp -p lib$(LIB)_nssa.a $(INSTALL_DIR)/lib/
	cp -p ./include_nssa/*.mod $(INSTALL_DIR)/include_nssa/
	@echo
	@echo

header: FORCE
	@echo
	@echo "================================"
	@echo ">>   Building common Physics    "
	@echo "================================"
	@echo
	@echo

mkdirs: FORCE
	@mkdir -p ./include_ns ./include_nssa ./include_ins ./include_ch ./include_mu ./include_slr ./include_slr_ins_v04
	@mkdir -p ./build_ns ./build_nssa ./build_ins ./build_ch ./build_mu ./build_slr ./build_slr_ins_v04

clean: FORCE
	$(RM) lib$(LIB)_ns.a lib$(LIB)_nssa.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_mu.a lib$(LIB)_slr.a lib$(LIB)_slr_ins_v04.a
	$(RM) ./build_ns/*.o ./include_ns/*.mod
	$(RM) ./build_nssa/*.o ./include_nssa/*.mod
	$(RM) ./build_ins/*.o ./include_ins/*.mod
	$(RM) ./build_ch/*.o ./include_ch/*.mod
	$(RM) ./build_mu/*.o ./include_mu/*.mod
	$(RM) ./build_slr_ins_v04/*.o ./include_slr_ins_v04/*.mod
	$(RM) ./*.mod

#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE: