#####################################################################
## HORSES - a High-Order Spectral Element Solver
##
##	Main code Makefile
##
##
#####################################################################

########################################
# 1.- User defined parameters
########################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIN?=ranlib
MOD_FLAG?=-J
PLATFORM?=LINUX
COMPILER?=gfortran
MODE?=RELEASE
COMM?=PARALLEL
MAKE=make
AR=ar
RANLIB=ranlib
PROG=./bin/horses3d.ns
PROG2=./bin/horses3d.nssa

OBJS_HORSES=	SpatialDiscretization \
		main  
OBJS_HORSESNSSA=	SpatialDiscretization \
		main  

OBJS_HORSES_BUILD=$(foreach obj,$(OBJS_HORSES),./build/$(obj).o)
OBJS_HORSESNSSA_BUILD=$(foreach obj,$(OBJS_HORSESNSSA),./build/$(obj).o)


##########################################
##	Include dependencies 		##
##########################################
-include ./make.deps

./build/main.o: FORCE


##########################################################
##		COMPILATION RULES 			##
##########################################################

.DEFAULT_GOAL := $(PROG)

./build/%.o: ./%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -I./include $(INCLUDE) $(INCLUDE)_ns $(MACROS) $(EXTLIB_INC) -c $< -o $@ $(MOD_FLAG) ./include

./build/%.o: ./%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -I./include $(INCLUDE) $(INCLUDE)_ns $(MACROS) $(EXTLIB_INC) -c $< -o $@ $(MOD_FLAG) ./include

./build_nssa/%.o: ./**/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DSPALARTALMARAS -I./include $(INCLUDE) $(INCLUDE)_nssa $(MACROS) $(EXTLIB_INC) -c $< -o $@ $(MOD_FLAG) ./include_nssa

./build_nssa/%.o: ./**/%.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DSPALARTALMARAS -I./include $(INCLUDE) $(INCLUDE)_nssa $(MACROS) $(EXTLIB_INC) -c $< -o $@ $(MOD_FLAG) ./include_nssa

$(PROG): header mkdirs $(OBJS_HORSES_BUILD)
	@echo
	@echo
	@echo "---------------------------------"
	@echo "| Linking " $@ "|"
	@echo "---------------------------------"
	$(FC) $(FFLAGS) $(INCLUDE) $(INCLUDE)_ns -o $(PROG) $(OBJS_HORSES_BUILD) $(HORSES_LIBS_EXEC) $(LIBS)

$(PROG): header mkdirs $(OBJS_HORSESNSSA_BUILD)
	@echo
	@echo
	@echo "---------------------------------"
	@echo "| Linking " $@ "|"
	@echo "---------------------------------"
	$(FC) $(FFLAGS) $(INCLUDE) $(INCLUDE)_nssa -o $(PROG2) $(OBJS_HORSESNSSA_BUILD) $(HORSES_LIBS_EXEC) $(LIBS)

install: $(PROG)
	cp -p $(PROG) $(INSTALL_DIR)/bin
	@echo
	@echo

install: $(PROG2)
	cp -p $(PROG2) $(INSTALL_DIR)/bin
	@echo
	@echo

##########################################
##	Extra auxiliar rules 		##
##########################################

clean: FORCE
	rm -f ./include/*.mod ./lib/*.a ./build/*.o
	rm -f $(PROG)
	rm -f ./include_nssa/*.mod ./lib_nssa/*.a ./build_nssa/*.o
	rm -f $(PROG2)
mkdirs: FORCE
	@mkdir -p ./build ./build_nssa
	@mkdir -p ./include ./include_nssa
	@mkdir -p ./bin ./build_nssa
	@mkdir -p ./lib ./lib_nssa

header: FORCE
	@echo
	@echo "============================================="
	@echo ">>   Building HORSES3D Navier-Stokes Solver  "
	@echo "============================================="
	@echo
	@echo

## Dummy target to make sure a rule is executed
FORCE:
