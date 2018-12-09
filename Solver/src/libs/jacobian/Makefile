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
			    ./build_ns/StaticCondensedMatrixClass.o \
			    ./build_ns/MatrixClass.o \
				./build_ns/NumericalJacobian.o

LIBOBJS_iNS = 	./build_ins/ColorsClass.o \
				./build_ins/Jacobian.o \
					./build_ins/GenericMatrixClass.o \
			    ./build_ins/LinkedListMatrixClass.o \
			    ./build_ins/CSRMatrixClass.o \
			    ./build_ins/PETScMatrixClass.o \
			    ./build_ins/DenseBlockDiagonalMatrixClass.o \
			    ./build_ins/SparseBlockDiagonalMatrixClass.o \
			    ./build_ins/StaticCondensedMatrixClass.o \
			    ./build_ins/MatrixClass.o \
				./build_ins/NumericalJacobian.o



LIBOBJS_CH = 	./build_ch/ColorsClass.o \
				  ./build_ch/Jacobian.o \
				  ./build_ch/GenericMatrixClass.o \
			    ./build_ch/LinkedListMatrixClass.o \
			    ./build_ch/CSRMatrixClass.o \
			    ./build_ch/PETScMatrixClass.o \
			    ./build_ch/DenseBlockDiagonalMatrixClass.o \
			    ./build_ch/SparseBlockDiagonalMatrixClass.o \
			    ./build_ch/StaticCondensedMatrixClass.o \
			    ./build_ch/MatrixClass.o \
				./build_ch/NumericalJacobian.o

LIBOBJS_NSCH = 	./build_nsch/ColorsClass.o \
				./build_nsch/Jacobian.o \
					./build_nsch/GenericMatrixClass.o \
			    ./build_nsch/LinkedListMatrixClass.o \
			    ./build_nsch/CSRMatrixClass.o \
			    ./build_nsch/PETScMatrixClass.o \
			    ./build_nsch/DenseBlockDiagonalMatrixClass.o \
			    ./build_nsch/SparseBlockDiagonalMatrixClass.o \
			    ./build_nsch/StaticCondensedMatrixClass.o \
			    ./build_nsch/MatrixClass.o \
				./build_nsch/NumericalJacobian.o

LIBOBJS_iNSCH =	./build_insch/ColorsClass.o \
				./build_insch/Jacobian.o \
					./build_insch/GenericMatrixClass.o \
			    ./build_insch/LinkedListMatrixClass.o \
			    ./build_insch/CSRMatrixClass.o \
			    ./build_insch/PETScMatrixClass.o \
			    ./build_insch/DenseBlockDiagonalMatrixClass.o \
			    ./build_insch/SparseBlockDiagonalMatrixClass.o \
			    ./build_insch/StaticCondensedMatrixClass.o \
			    ./build_insch/MatrixClass.o \
				./build_insch/NumericalJacobian.o



LIB = jacobian

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
	@echo ">>   Building Jacobian    "
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
