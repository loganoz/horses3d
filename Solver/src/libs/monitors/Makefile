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
		./build_ns/Monitors.o

LIBOBJS_iNS = 	./build_ins/MonitorDefinitions.o \
		./build_ins/Probe.o \
		./build_ins/StatisticsMonitor.o \
		./build_ins/SurfaceIntegrals.o \
		./build_ins/VolumeIntegrals.o \
		./build_ins/ResidualsMonitor.o \
		./build_ins/SurfaceMonitor.o \
		./build_ins/VolumeMonitor.o \
		./build_ins/Monitors.o

LIBOBJS_CH = 	./build_ch/MonitorDefinitions.o \
		./build_ch/Probe.o\
		./build_ch/StatisticsMonitor.o \
		./build_ch/SurfaceIntegrals.o \
		./build_ch/VolumeIntegrals.o \
		./build_ch/ResidualsMonitor.o \
		./build_ch/SurfaceMonitor.o \
		./build_ch/VolumeMonitor.o \
		./build_ch/Monitors.o

LIBOBJS_NSCH = 	./build_nsch/MonitorDefinitions.o \
		./build_nsch/Probe.o\
		./build_nsch/StatisticsMonitor.o \
		./build_nsch/SurfaceIntegrals.o \
		./build_nsch/VolumeIntegrals.o \
		./build_nsch/ResidualsMonitor.o \
		./build_nsch/SurfaceMonitor.o \
		./build_nsch/VolumeMonitor.o \
		./build_nsch/Monitors.o

LIBOBJS_iNSCH = ./build_insch/MonitorDefinitions.o \
		./build_insch/Probe.o\
		./build_insch/StatisticsMonitor.o \
		./build_insch/SurfaceIntegrals.o \
		./build_insch/VolumeIntegrals.o \
		./build_insch/ResidualsMonitor.o \
		./build_insch/SurfaceMonitor.o \
		./build_insch/VolumeMonitor.o \
		./build_insch/Monitors.o


LIB = monitors

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
	@echo ">>   Building Monitors    "
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

