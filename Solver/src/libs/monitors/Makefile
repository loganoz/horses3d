#################################################
#	Monitors Makefile 		#
#################################################
INSTALL_DIR?=./
INCLUDE?=-I../../../include
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J


LIBOBJS_NS = 	./build_ns/MonitorDefinitions.o \
		./build_ns/Probe.o \
		./build_ns/StatisticsMonitor.o \
		./build_ns/SurfaceIntegrals.o \
		./build_ns/VolumeIntegrals.o \
		./build_ns/ResidualsMonitor.o \
		./build_ns/SurfaceMonitor.o \
		./build_ns/VolumeMonitor.o \
		./build_ns/LoadBalancingMonitor.o \
		./build_ns/Monitors.o \
		./build_ns/FWHDefinitions.o \
		./build_ns/FWHObseverClass.o \
		./build_ns/FWHGeneralClass.o

LIBOBJS_NSSA = 	./build_nssa/MonitorDefinitions.o \
		./build_nssa/Probe.o \
		./build_nssa/StatisticsMonitor.o \
		./build_nssa/SurfaceIntegrals.o \
		./build_nssa/VolumeIntegrals.o \
		./build_nssa/ResidualsMonitor.o \
		./build_nssa/SurfaceMonitor.o \
		./build_nssa/VolumeMonitor.o \
		./build_nssa/LoadBalancingMonitor.o \
		./build_nssa/Monitors.o

LIBOBJS_iNS = 	./build_ins/MonitorDefinitions.o \
		./build_ins/Probe.o \
		./build_ins/StatisticsMonitor.o \
		./build_ins/SurfaceIntegrals.o \
		./build_ins/VolumeIntegrals.o \
		./build_ins/ResidualsMonitor.o \
		./build_ins/SurfaceMonitor.o \
		./build_ins/VolumeMonitor.o \
		./build_ins/LoadBalancingMonitor.o \
		./build_ins/Monitors.o

LIBOBJS_CH = 	./build_ch/MonitorDefinitions.o \
		./build_ch/Probe.o\
		./build_ch/StatisticsMonitor.o \
		./build_ch/SurfaceIntegrals.o \
		./build_ch/VolumeIntegrals.o \
		./build_ch/ResidualsMonitor.o \
		./build_ch/SurfaceMonitor.o \
		./build_ch/VolumeMonitor.o \
		./build_ch/LoadBalancingMonitor.o \
		./build_ch/Monitors.o

LIBOBJS_MU = 	./build_mu/MonitorDefinitions.o \
		./build_mu/Probe.o\
		./build_mu/StatisticsMonitor.o \
		./build_mu/SurfaceIntegrals.o \
		./build_mu/VolumeIntegrals.o \
		./build_mu/ResidualsMonitor.o \
		./build_mu/SurfaceMonitor.o \
		./build_mu/VolumeMonitor.o \
		./build_mu/LoadBalancingMonitor.o \
		./build_mu/Monitors.o

LIB = monitors

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:=all

all: lib$(LIB)_ns.a lib$(LIB)_nssa.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_mu.a

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

install: all
	cp -p lib$(LIB)_ns.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_nssa.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_ins.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_ch.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_mu.a $(INSTALL_DIR)/lib/
	cp -p ./include_ns/*.mod $(INSTALL_DIR)/include_ns/
	cp -p ./include_nssa/*.mod $(INSTALL_DIR)/include_nssa/
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

install-nssa: lib$(LIB)_nssa.a
	cp -p lib$(LIB)_nssa.a $(INSTALL_DIR)/lib/
	cp -p ./include_nssa/*.mod $(INSTALL_DIR)/include_nssa/
	@echo
	@echo

header: FORCE
	@echo
	@echo "================================"
	@echo ">>   Building Monitors    "
	@echo "================================"
	@echo
	@echo

mkdirs: FORCE
	@mkdir -p ./include_ns ./include_nssa ./include_ins ./include_ch ./include_mu
	@mkdir -p ./build_ns ./build_nssa ./build_ins ./build_ch ./build_mu

clean: FORCE
	$(RM) lib$(LIB)_ns.a lib$(LIB)_nssa.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_mu.a
	$(RM) ./build_ns/*.o ./include_ns/*.mod
	$(RM) ./build_nssa/*.o ./include_nssa/*.mod
	$(RM) ./build_ins/*.o ./include_ins/*.mod
	$(RM) ./build_ch/*.o ./include_ch/*.mod
	$(RM) ./build_mu/*.o ./include_mu/*.mod

#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE: