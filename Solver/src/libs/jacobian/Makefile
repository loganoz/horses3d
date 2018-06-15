#################################################
#	Jacobian Makefile 		#
#################################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J




LIBOBJS_NS = 	./build_ns/ColorsClass.o \
				./build_ns/Jacobian.o \
					./build_ns/GenericMatrixClass.o \
			    ./build_ns/LinkedListMatrixClass.o \
			    ./build_ns/CSRMatrixClass.o \
			    ./build_ns/PETScMatrixClass.o \
			    ./build_ns/DenseBlockDiagonalMatrixClass.o \
			    ./build_ns/SparseBlockDiagonalMatrixClass.o \
			    ./build_ns/MatrixClass.o \
				./build_ns/NumericalJacobian.o

LIBOBJS_CH = 	./build_ch/ColorsClass.o \
				  ./build_ch/Jacobian.o \
				  ./build_ch/GenericMatrixClass.o \
			    ./build_ch/LinkedListMatrixClass.o \
			    ./build_ch/CSRMatrixClass.o \
			    ./build_ch/PETScMatrixClass.o \
			    ./build_ch/DenseBlockDiagonalMatrixClass.o \
			    ./build_ch/SparseBlockDiagonalMatrixClass.o \
			    ./build_ch/MatrixClass.o \
				./build_ch/NumericalJacobian.o

LIBOBJS_MU = 	./build_mu/ColorsClass.o \
				./build_mu/Jacobian.o \
					./build_mu/GenericMatrixClass.o \
			    ./build_mu/LinkedListMatrixClass.o \
			    ./build_mu/CSRMatrixClass.o \
			    ./build_mu/PETScMatrixClass.o \
			    ./build_mu/DenseBlockDiagonalMatrixClass.o \
			    ./build_mu/SparseBlockDiagonalMatrixClass.o \
			    ./build_mu/MatrixClass.o \
				./build_mu/NumericalJacobian.o


LIB = jacobian

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:=all

all: lib$(LIB)_ns.a lib$(LIB)_ch.a lib$(LIB)_mu.a

./build_ns/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -I./include_ns $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ns -c $< -o $@ $(MOD_FLAG) ./include_ns

./build_ch/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DCAHNHILLIARD -I./include_ch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ch -c $< -o $@ $(MOD_FLAG) ./include_ch

./build_mu/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -DCAHNHILLIARD -I./include_mu $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_mu -c $< -o $@ $(MOD_FLAG) ./include_mu

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
	cp -p lib$(LIB)_ch.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_mu.a $(INSTALL_DIR)/lib/
	cp -p ./include_ns/*.mod $(INSTALL_DIR)/include_ns/
	cp -p ./include_ch/*.mod $(INSTALL_DIR)/include_ch/
	cp -p ./include_mu/*.mod $(INSTALL_DIR)/include_mu/
	@echo
	@echo

header: FORCE
	@echo
	@echo "================================"
	@echo ">>   Building Jacobian    "
	@echo "================================"
	@echo
	@echo

mkdirs: FORCE
	@mkdir -p ./include_ns ./include_ch ./include_mu
	@mkdir -p ./build_ns ./build_ch ./build_mu

clean: FORCE
	$(RM) lib$(LIB)_ns.a  lib$(LIB)_ch.a
	$(RM) ./build_ns/*.o ./include_ns/*.mod
	$(RM) ./build_ch/*.o ./include_ch/*.mod
	$(RM) ./build_mu/*.o ./include_mu/*.mod

#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE:
