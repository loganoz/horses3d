#################################################
#	Spectral Makefile 		#
#################################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J

LIBOBJS_NS =    ./build_ns/DGSEMClass.o \
		./build_ns/DGIntegrals.o \
		./build_ns/GradientsStabilization.o \
		./build_ns/HyperbolicDiscretizationClass.o \
		./build_ns/HyperbolicStandard.o \
		./build_ns/HyperbolicSplitForm.o \
		./build_ns/HyperbolicDiscretizations.o \
		./build_ns/EllipticDiscretizationClass.o \
		./build_ns/EllipticBR1.o \
		./build_ns/EllipticBR2.o \
		./build_ns/EllipticIP.o \
		./build_ns/EllipticDiscretizations.o \
		./build_ns/SVV.o 

LIBOBJS_iNS =    ./build_ins/DGSEMClass.o \
		./build_ins/DGIntegrals.o \
		./build_ins/GradientsStabilization.o \
		./build_ins/HyperbolicDiscretizationClass.o \
		./build_ins/HyperbolicStandard.o \
		./build_ins/HyperbolicSplitForm.o \
		./build_ins/HyperbolicDiscretizations.o \
		./build_ins/EllipticDiscretizationClass.o \
		./build_ins/EllipticBR1.o \
		./build_ins/EllipticBR2.o \
		./build_ins/EllipticIP.o \
		./build_ins/EllipticDiscretizations.o \
		./build_ins/SVV.o 

LIBOBJS_CH =    ./build_ch/DGSEMClass.o \
		./build_ch/DGIntegrals.o \
		./build_ch/GradientsStabilization.o \
		./build_ch/HyperbolicDiscretizationClass.o \
		./build_ch/HyperbolicStandard.o \
		./build_ch/HyperbolicSplitForm.o \
		./build_ch/HyperbolicDiscretizations.o \
		./build_ch/EllipticDiscretizationClass.o \
		./build_ch/EllipticBR1.o \
		./build_ch/EllipticBR2.o \
		./build_ch/EllipticIP.o \
		./build_ch/EllipticDiscretizations.o \
		./build_ch/SVV.o 

LIBOBJS_NSCH =  ./build_nsch/DGSEMClass.o \
		./build_nsch/DGIntegrals.o \
		./build_nsch/GradientsStabilization.o \
		./build_nsch/HyperbolicDiscretizationClass.o \
		./build_nsch/HyperbolicStandard.o \
		./build_nsch/HyperbolicSplitForm.o \
		./build_nsch/HyperbolicDiscretizations.o \
		./build_nsch/EllipticDiscretizationClass.o \
		./build_nsch/EllipticBR1.o \
		./build_nsch/EllipticBR2.o \
		./build_nsch/EllipticIP.o \
		./build_nsch/EllipticDiscretizations.o \
		./build_nsch/SVV.o 

LIBOBJS_iNSCH = ./build_insch/DGSEMClass.o \
		./build_insch/DGIntegrals.o \
		./build_insch/GradientsStabilization.o \
		./build_insch/HyperbolicDiscretizationClass.o \
		./build_insch/HyperbolicStandard.o \
		./build_insch/HyperbolicSplitForm.o \
		./build_insch/HyperbolicDiscretizations.o \
		./build_insch/EllipticDiscretizationClass.o \
		./build_insch/EllipticBR1.o \
		./build_insch/EllipticBR2.o \
		./build_insch/EllipticIP.o \
		./build_insch/EllipticDiscretizations.o \
		./build_insch/SVV.o 

LIB = discretization

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:=all

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
	@echo ">>   Building Discretization    "
	@echo "================================"
	@echo
	@echo

mkdirs: FORCE
	@mkdir -p ./include_ns ./include_ins ./include_ch ./include_nsch ./include_insch
	@mkdir -p ./build_ns ./build_ins ./build_ch ./build_nsch ./build_insch

clean: FORCE
	$(RM) lib$(LIB)_ns.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_nsch.a lib$(LIB)_insch.a
	$(RM) ./build_ns/*.o ./include_ns/*.mod
	$(RM) ./build_ins/*.o ./include_ins/*.mod
	$(RM) ./build_ch/*.o ./include_ch/*.mod
	$(RM) ./build_nsch/*.o ./include_nsch/*.mod
	$(RM) ./build_insch/*.o ./include_insch/*.mod

#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE:
