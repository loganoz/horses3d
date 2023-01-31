#################################################
#	Foundation Makefile 		#
#################################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J 

LIBOBJS = ./build/FluidData_NSSA.o \
		  ./build/PhysicsStorage_NSSA.o \
		  ./build/VariableConversion_NSSA.o \
		  ./build/ArtificialViscosity.o \
		  ./build/Turbulence.o  \
		  ./build/Physics_NSSA.o \
		  ./build/RiemannSolvers_NSSA.o \
		  ./build/ManufacturedSolutions.o\

          
LIB = physicsnssa

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:=  lib$(LIB).a


./build/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -DSPALARTALMARAS $(MACROS) $(EXTLIB_INC) -I./include $(INCLUDE) -c $< -o $@ $(MOD_FLAG) ./include

lib$(LIB).a: header mkdirs $(LIBOBJS)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS)
	@$(RANLIB) $@

install : lib$(LIB).a
	cp -p lib$(LIB).a $(INSTALL_DIR)/lib/
	cp -p ./include/*.mod $(INSTALL_DIR)/include/
	@echo
	@echo

header: FORCE
	@echo
	@echo "============================================================"
	@echo ">>   Building Navier-Stokes/Spallart-Almaras physics        "
	@echo "============================================================"
	@echo
	@echo
mkdirs: FORCE
	@mkdir -p ./include
	@mkdir -p ./build

clean: FORCE
	$(RM) lib$(LIB).a 
	$(RM) ./build/*.o ./include/*.mod

#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE: