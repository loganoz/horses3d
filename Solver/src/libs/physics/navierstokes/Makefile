#################################################
#	Foundation Makefile 		#
#################################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J




LIBOBJS = ./build/FluidData_NS.o \
		  ./build/PhysicsStorage_NS.o \
		  ./build/VariableConversion_NS.o \
		  ./build/ArtificialViscosity.o \
		  ./build/LESModels.o \
		  ./build/Physics_NS.o \
		  ./build/RiemannSolvers_NS.o \
		  ./build/ManufacturedSolutions.o \

LIBOBJS_NSSA = 	./build/FluidData_NS.o \
		  ./build/PhysicsStorage_NS.o \
		  ./build/VariableConversion_NS.o \
		  ./build/ArtificialViscosity.o \
		  ./build/LESModels.o \
		  ./build/Physics_NS.o \
		  ./build/RiemannSolvers_NS.o \
		  ./build/ManufacturedSolutions.o \
		  ./build/Turbulence.o  
          
LIB = physicsns

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL := lib$(LIB).a

all: lib$(LIB).a lib$(LIB)_nssa.a 


./build/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(MACROS) $(EXTLIB_INC) -I./include $(INCLUDE) -c $< -o $@ $(MOD_FLAG) ./include

./build_nssa/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DSPALARTALMARAS $(MACROS) $(EXTLIB_INC) -I./include_nssa $(INCLUDE) $(INCLUDE)_nssa -c $< -o $@ $(MOD_FLAG) ./include_nssa

lib$(LIB).a: header mkdirs $(LIBOBJS)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS)
	@$(RANLIB) $@

lib$(LIB)_nssa.a: header mkdirs $(LIBOBJS_NSSA)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_NSSA)
	@$(RANLIB) $@

install: lib$(LIB).a
	cp -p lib$(LIB).a $(INSTALL_DIR)/lib/
	cp -p ./include/*.mod $(INSTALL_DIR)/include/
	@echo
	@echo

install: lib$(LIB)_nssa.a
	cp -p lib$(LIB)_nssa.a $(INSTALL_DIR)/lib/
	cp -p ./include_nssa/*.mod $(INSTALL_DIR)/include_nssa/
	@echo
	@echo

header: FORCE
	@echo
	@echo "==========================================="
	@echo ">>   Building Navier-Stokes physics        "
	@echo "==========================================="
	@echo
	@echo
mkdirs: FORCE
	@mkdir -p ./include ./include_nssa
	@mkdir -p ./build ./build_nssa

clean: FORCE
	$(RM) lib$(LIB).a lib$(LIB)_nssa.a
	$(RM) ./build/*.o ./include/*.mod
	$(RM) ./build/_nssa/*.o ./include_nssa/*.mod

#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE:
