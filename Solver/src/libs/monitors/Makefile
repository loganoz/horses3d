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
		./build_ns/Probe.o\
		./build_ns/StatisticsMonitor.o \
		./build_ns/SurfaceIntegrals.o \
		./build_ns/VolumeIntegrals.o \
		./build_ns/ResidualsMonitor.o \
		./build_ns/SurfaceMonitor.o \
		./build_ns/VolumeMonitor.o \
		./build_ns/Monitors.o

LIBOBJS_CH = 	./build_ch/MonitorDefinitions.o \
		./build_ch/Probe.o\
		./build_ch/StatisticsMonitor.o \
		./build_ch/SurfaceIntegrals.o \
		./build_ch/VolumeIntegrals.o \
		./build_ch/ResidualsMonitor.o \
		./build_ch/SurfaceMonitor.o \
		./build_ch/VolumeMonitor.o \
		./build_ch/Monitors.o

LIB = monitors

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:=all

all: lib$(LIB)_ns.a

./build_ns/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -I./include_ns $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ns -c $< -o $@ $(MOD_FLAG) ./include_ns

./build_ch/%.o: %.f90
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
	@echo ">>   Building Monitors    "
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

