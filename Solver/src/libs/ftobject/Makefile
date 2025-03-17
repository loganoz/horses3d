#################################################
#	FTObjectLibrary Makefile 		#
#################################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J




LIBOBJS = 	./build/FTOLConstants.o \
	 	./build/FTObjectClass.o \
		./build/FTDataClass.o \
		./build/FTObjectArrayClass.o \
		./build/FTLinkedListClass.o \
	 	./build/FTMultiIndexTable.o \
		./build/FTDictionaryClass.o \
		./build/FTValueClass.o \
		./build/FTValueDictionaryClass.o \
	 	./build/FTStackClass.o \
		./build/Hash.o \
		./build/FTSparseMatrixClass.o \
		./build/FTExceptionClass.o \
	 	./build/Assert.o \
		./build/Comparisons.o \
		./build/TestSuiteManagerClass.o \
	 	./build/FTObjectLibrary.o


LIB = ftobject


-include ./make.deps

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL := lib$(LIB).a

./build/%.o: %.f90
	@echo
	@echo $<
	$(FC) -fPIC $(FFLAGS) $(MACROS) $(EXTLIB_INC) -I ./include $(INCLUDE) -c $< -o $@ $(MOD_FLAG) ./include

lib$(LIB).a: header mkdirs $(LIBOBJS)
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS)
	@$(RANLIB) $@

install: lib$(LIB).a
	cp -p lib$(LIB).a $(INSTALL_DIR)/lib/
	cp -p ./include/*.mod $(INSTALL_DIR)/include/
	@echo
	@echo

header: FORCE
	@echo
	@echo "================================"
	@echo ">>   Building FTObjectLibrary   "
	@echo "================================"
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