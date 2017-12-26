#################################################
#	Mesh Makefile 		#
#################################################
INSTALL_DIR?=./
FFLAGS?=-cpp
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J

LIBOBJS = 	./build/ConnectivityClass.o \
		./build/MeshTypes.o \
		./build/StorageClass.o \
		./build/FacePatchClass.o \
		./build/TransfiniteMaps3D.o \
		./build/MappedGeometry.o \
		./build/NodeClass.o \
		./build/HexElementConnectivityDefinitions.o \
		./build/FaceClass.o \
		./build/HexElementClass.o \
		./build/ZoneClass.o \
		./build/HexMesh.o \
		./build/MeshPartitioning.o \
		./build/Read_HDF5Mesh_HOPR.o \
		./build/Read_SpecMesh.o \
		./build/ReadMeshFile.o

LIB = mesh

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL := lib$(LIB).a

./build/METISPartitioning.o: ./METISPartitioning.f90		# METIS calls require this, and only this flags... so it is separated
	@echo
	@echo $<
	$(FC) -cpp -O2 -ffree-line-length-0 -I./include $(MACROS) $(EXTLIB_INC) $(INCLUDE) -c $< -o $@ $(MOD_FLAG) ./include
	
./build/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -I./include $(MACROS) $(EXTLIB_INC) $(INCLUDE) -c $< -o $@ $(MOD_FLAG) ./include

lib$(LIB).a: header mkdirs $(LIBOBJS) ./build/METISPartitioning.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS) ./build/METISPartitioning.o
	@$(RANLIB) $@

install: lib$(LIB).a
	cp -p lib$(LIB).a $(INSTALL_DIR)/lib/
	cp -p ./include/*.mod $(INSTALL_DIR)/include/
	@echo
	@echo

header: FORCE
	@echo
	@echo "================================"
	@echo ">>   Building Mesh        "
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
