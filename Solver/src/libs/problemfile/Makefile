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

ifeq ($(PLATFORM),LINUX)
   LIBNAME=$(LIB).so
else ifeq ($(PLATFORM),MACOSX)
   LIBNAME=$(LIB).dylib
endif



############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL := $(LIBNAME)

./build/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) $(INCLUDE) -c $< -o $@ $(MOD_FLAG) ./include

./build/ProblemFile.o: ./ProblemFile.f90
	$(FC) -fPIC -shared $(FFLAGS) $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(MOD_FLAG)./include -c $< -o $@

lib$(LIBNAME): header mkdirs ./build/ProblemFile.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	$(FC) $(DYNLIB_FLAG) $(FFLAGS) $(MACROS) $(INCLUDE) ./build/ProblemFile.o $(LIBS)

install: lib$(LIBNAME)
	cp -p lib$(LIBNAME) $(INSTALL_DIR)/lib/
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
	@mkdir -p ./build

clean: FORCE
	$(RM) $(LIBNAME)
	$(RM) ./build/*.o ./include/*.mod

#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE:
