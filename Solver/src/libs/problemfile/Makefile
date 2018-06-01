#################################################
#	Foundation Makefile 		#
#################################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J
PLATFORM?=LINUX


LIB = problemfile

ifeq ($(PLATFORM),MACOSX)
   DYNLIB_FLAG_NS=-fPIC -shared -install_name @rpath/libproblemfile_ns.dylib -o ./libproblemfile_ns.dylib
   DYNLIB_FLAG_CH=-fPIC -shared -install_name @rpath/libproblemfile_ch.dylib -o ./libproblemfile_ch.dylib
   DYNLIB_FLAG_MU=-fPIC -shared -install_name @rpath/libproblemfile_mu.dylib -o ./libproblemfile_mu.dylib
   LIBSUFFIX=.dylib
else ifeq ($(PLATFORM),LINUX)
   DYNLIB_FLAG_NS=-fPIC -shared -o ./libproblemfile_ns.so
   DYNLIB_FLAG_CH=-fPIC -shared -o ./libproblemfile_ch.so
   DYNLIB_FLAG_MU=-fPIC -shared -o ./libproblemfile_mu.so
   LIBSUFFIX=.so
endif




############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:=all

all: lib$(LIB)_ns$(LIBSUFFIX) lib$(LIB)_ch$(LIBSUFFIX) lib$(LIB)_mu$(LIBSUFFIX)

./build_ns/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(INCLUDE) -c $< -o $@ $(MOD_FLAG) ./include

./build_ch/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(INCLUDE) -c $< -o $@ $(MOD_FLAG) ./include

./build_mu/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(INCLUDE) -c $< -o $@ $(MOD_FLAG) ./include

./build_ns/ProblemFile.o: ./ProblemFile.f90
	$(FC) -fPIC -shared -DNAVIERSTOKES $(FFLAGS) $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ns $(MOD_FLAG)./include -c $< -o $@

./build_ch/ProblemFile.o: ./ProblemFile.f90
	$(FC) -fPIC -shared -DCAHNHILLIARD $(FFLAGS) $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ch $(MOD_FLAG)./include -c $< -o $@

./build_mu/ProblemFile.o: ./ProblemFile.f90
	$(FC) -fPIC -shared -DNAVIERSTOKES -DCAHNHILLIARD $(FFLAGS) $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_mu $(MOD_FLAG)./include -c $< -o $@

lib$(LIB)_ns$(LIBSUFFIX): header mkdirs ./build_ns/ProblemFile.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	$(FC) $(DYNLIB_FLAG_NS) $(FFLAGS) $(MACROS) $(INCLUDE)_ns ./build_ns/ProblemFile.o $(LIBS)

lib$(LIB)_ch$(LIBSUFFIX): header mkdirs ./build_ch/ProblemFile.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	$(FC) $(DYNLIB_FLAG_CH) $(FFLAGS) $(MACROS) $(INCLUDE)_ch ./build_ch/ProblemFile.o $(LIBS)

lib$(LIB)_mu$(LIBSUFFIX): header mkdirs ./build_mu/ProblemFile.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	$(FC) $(DYNLIB_FLAG_MU) $(FFLAGS) $(MACROS) $(INCLUDE)_mu ./build_mu/ProblemFile.o $(LIBS)

install: all
	cp -p lib$(LIB)_ns$(LIBSUFFIX) $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_ch$(LIBSUFFIX) $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_mu$(LIBSUFFIX) $(INSTALL_DIR)/lib/
	@echo
	@echo

header: FORCE
	@echo
	@echo "================================"
	@echo ">>   Building Problem File      "
	@echo "================================"
	@echo
	@echo

mkdirs: FORCE
	@mkdir -p ./include
	@mkdir -p ./build_ns ./build_ch ./build_mu

clean: FORCE
	$(RM) $(LIBNAME)
	$(RM) ./build_ns/*.o ./build_ch/*.o ./build_mu/*.o ./include/*.mod

#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE:
