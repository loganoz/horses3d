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
		./build_ns/BoundaryConditions.o 

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
		./build_ch/BoundaryConditions.o 

LIBOBJS_NSCH = 	./build_nsch/FluidData.o \
		./build_nsch/PhysicsStorage.o \
		./build_nsch/VariableConversion.o \
		./build_nsch/Physics.o \
		./build_nsch/GenericBoundaryConditionClass.o \
		./build_nsch/InflowBC.o \
		./build_nsch/OutflowBC.o \
		./build_nsch/NoSlipWallBC.o \
		./build_nsch/FreeSlipWallBC.o \
		./build_nsch/PeriodicBC.o \
		./build_nsch/UserDefinedBC.o \
		./build_nsch/BoundaryConditions.o 

LIBOBJS_iNSCH =	./build_insch/FluidData.o \
		./build_insch/PhysicsStorage.o \
		./build_insch/VariableConversion.o \
		./build_insch/Physics.o \
		./build_insch/GenericBoundaryConditionClass.o \
		./build_insch/InflowBC.o \
		./build_insch/OutflowBC.o \
		./build_insch/NoSlipWallBC.o \
		./build_insch/FreeSlipWallBC.o \
		./build_insch/PeriodicBC.o \
		./build_insch/UserDefinedBC.o \
		./build_insch/BoundaryConditions.o 

LIB = physics

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL := lib$(LIB).a

all: lib$(LIB)_ns.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_nsch.a lib$(LIB)_insch.a

./build_ns/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -I./include_ns $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ns -c $< -o $@ $(MOD_FLAG) ./include_ns

./build_ins/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DINCNS -I./include_ins $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ins -c $< -o $@ $(MOD_FLAG) ./include_ins

./build_ch/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DCAHNHILLIARD -I./include_ch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ch -c $< -o $@ $(MOD_FLAG) ./include_ch

./build_nsch/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -DCAHNHILLIARD -I./include_nsch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nsch -c $< -o $@ $(MOD_FLAG) ./include_nsch

./build_insch/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DINCNS -DCAHNHILLIARD -I./include_insch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_insch -c $< -o $@ $(MOD_FLAG) ./include_insch

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

lib$(LIB)_nsch.a: header mkdirs $(LIBOBJS_NSCH)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_NSCH) 
	@$(RANLIB) $@

lib$(LIB)_insch.a: header mkdirs $(LIBOBJS_iNSCH)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_iNSCH) 
	@$(RANLIB) $@

install: all
	cp -p lib$(LIB)_ns.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_ins.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_ch.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_nsch.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_insch.a $(INSTALL_DIR)/lib/
	cp -p ./include_ns/*.mod $(INSTALL_DIR)/include_ns/
	cp -p ./include_ins/*.mod $(INSTALL_DIR)/include_ins/
	cp -p ./include_ch/*.mod $(INSTALL_DIR)/include_ch/
	cp -p ./include_nsch/*.mod $(INSTALL_DIR)/include_nsch/
	cp -p ./include_insch/*.mod $(INSTALL_DIR)/include_insch/
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
	@mkdir -p ./include_ns ./include_ins ./include_ch ./include_nsch ./include_insch
	@mkdir -p ./build_ns ./build_ins ./build_ch ./build_nsch ./build_insch

clean: FORCE
	$(RM) lib$(LIB)_ns.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_nsch.a
	$(RM) ./build_ns/*.o ./include_ns/*.mod
	$(RM) ./build_ins/*.o ./include_ins/*.mod
	$(RM) ./build_ch/*.o ./include_ch/*.mod
	$(RM) ./build_nsch/*.o ./include_nsch/*.mod
	$(RM) ./build_insch/*.o ./include_insch/*.mod

#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE:
