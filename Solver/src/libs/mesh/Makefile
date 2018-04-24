#################################################
#	Mesh Makefile 		#
#################################################
INSTALL_DIR?=./
FFLAGS?=-cpp -ffree-line-length-0
AR?=ar
ARFLAGS?=cr
RANLIB?=ranlib
MOD_FLAG?=-J

LIBOBJS_NS = 	./build_ns/ConnectivityClass.o \
		./build_ns/MeshTypes.o \
		./build_ns/WallDistance.o \
		./build_ns/StorageClass.o \
		./build_ns/FacePatchClass.o \
		./build_ns/TransfiniteMaps3D.o \
		./build_ns/MappedGeometry.o \
		./build_ns/NodeClass.o \
		./build_ns/HexElementConnectivityDefinitions.o \
		./build_ns/FaceClass.o \
		./build_ns/HexElementClass.o \
		./build_ns/ZoneClass.o \
		./build_ns/HexMesh.o \
		./build_ns/MeshPartitioning.o \
		./build_ns/Read_HDF5Mesh_HOPR.o \
		./build_ns/Read_SpecMesh.o \
		./build_ns/ReadMeshFile.o

LIBOBJS_CH = 	./build_ch/ConnectivityClass.o \
		./build_ch/MeshTypes.o \
		./build_ch/WallDistance.o \
		./build_ch/StorageClass.o \
		./build_ch/FacePatchClass.o \
		./build_ch/TransfiniteMaps3D.o \
		./build_ch/MappedGeometry.o \
		./build_ch/NodeClass.o \
		./build_ch/HexElementConnectivityDefinitions.o \
		./build_ch/FaceClass.o \
		./build_ch/HexElementClass.o \
		./build_ch/ZoneClass.o \
		./build_ch/HexMesh.o \
		./build_ch/MeshPartitioning.o \
		./build_ch/Read_HDF5Mesh_HOPR.o \
		./build_ch/Read_SpecMesh.o \
		./build_ch/ReadMeshFile.o

LIBOBJS_NSCH =  ./build_nsch/ConnectivityClass.o \
		./build_nsch/MeshTypes.o \
		./build_nsch/WallDistance.o \
		./build_nsch/StorageClass.o \
		./build_nsch/FacePatchClass.o \
		./build_nsch/TransfiniteMaps3D.o \
		./build_nsch/MappedGeometry.o \
		./build_nsch/NodeClass.o \
		./build_nsch/HexElementConnectivityDefinitions.o \
		./build_nsch/FaceClass.o \
		./build_nsch/HexElementClass.o \
		./build_nsch/ZoneClass.o \
		./build_nsch/HexMesh.o \
		./build_nsch/MeshPartitioning.o \
		./build_nsch/Read_HDF5Mesh_HOPR.o \
		./build_nsch/Read_SpecMesh.o \
		./build_nsch/ReadMeshFile.o



LIB = mesh

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:=all

all: lib$(LIB)_ns.a lib$(LIB)_ch.a lib$(LIB)_nsch.a

./build_ns/METISPartitioning.o: ./METISPartitioning.f90		# METIS calls require this, and only this flags... so it is separated
	@echo
	@echo $<
	$(FC) -cpp -O2 -ffree-line-length-0 -DNAVIERSTOKES -I./include_ns $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ns -c $< -o $@ $(MOD_FLAG) ./include_ns

./build_ch/METISPartitioning.o: ./METISPartitioning.f90		# METIS calls require this, and only this flags... so it is separated
	@echo
	@echo $<
	$(FC) -cpp -O2 -ffree-line-length-0 -DCAHNHILLIARD -I./include_ch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ch -c $< -o $@ $(MOD_FLAG) ./include_ch

./build_nsch/METISPartitioning.o: ./METISPartitioning.f90		# METIS calls require this, and only this flags... so it is separated
	@echo
	@echo $<
	$(FC) -cpp -O2 -ffree-line-length-0 -DNAVIERSTOKES -DCAHNHILLIARD -I./include_nsch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nsch -c $< -o $@ $(MOD_FLAG) ./include_nsch

./build_ns/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -I./include_ns $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ns -c $< -o $@ $(MOD_FLAG) ./include_ns

./build_ch/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DCAHNHILLIARD -I./include_ch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ch -c $< -o $@ $(MOD_FLAG) ./include_ch

./build_nsch/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DNAVIERSTOKES -DCAHNHILLIARD -I./include_nsch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nsch -c $< -o $@ $(MOD_FLAG) ./include_nsch

lib$(LIB)_ns.a: header mkdirs $(LIBOBJS_NS) ./build_ns/METISPartitioning.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_NS) ./build_ns/METISPartitioning.o
	@$(RANLIB) $@

lib$(LIB)_ch.a: header mkdirs $(LIBOBJS_CH) ./build_ch/METISPartitioning.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_CH) ./build_ch/METISPartitioning.o
	@$(RANLIB) $@

lib$(LIB)_nsch.a: header mkdirs $(LIBOBJS_NSCH) ./build_nsch/METISPartitioning.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_NSCH) ./build_nsch/METISPartitioning.o
	@$(RANLIB) $@

install: all
	cp -p lib$(LIB)_ns.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_ch.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_nsch.a $(INSTALL_DIR)/lib/
	cp -p ./include_ns/*.mod $(INSTALL_DIR)/include_ns/
	cp -p ./include_ch/*.mod $(INSTALL_DIR)/include_ch/
	cp -p ./include_nsch/*.mod $(INSTALL_DIR)/include_nsch/
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
	@mkdir -p ./include_ns ./include_ch ./include_nsch
	@mkdir -p ./build_ns ./build_ch ./build_nsch

clean: FORCE
	$(RM) lib$(LIB)_ns.a  lib$(LIB)_ch.a
	$(RM) ./build_ns/*.o ./include_ns/*.mod 
	$(RM) ./build_ch/*.o ./include_ch/*.mod
	$(RM) ./build_nsch/*.o ./include_nsch/*.mod

#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE:
