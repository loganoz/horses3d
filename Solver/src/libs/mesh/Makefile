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
		./build_ns/Read_GMSH.o \
		./build_ns/ReadMeshFile.o

LIBOBJS_iNS = 	./build_ins/ConnectivityClass.o \
		./build_ins/MeshTypes.o \
		./build_ins/WallDistance.o \
		./build_ins/StorageClass.o \
		./build_ins/FacePatchClass.o \
		./build_ins/TransfiniteMaps3D.o \
		./build_ins/MappedGeometry.o \
		./build_ins/NodeClass.o \
		./build_ins/HexElementConnectivityDefinitions.o \
		./build_ins/FaceClass.o \
		./build_ins/HexElementClass.o \
		./build_ins/ZoneClass.o \
		./build_ins/HexMesh.o \
		./build_ins/MeshPartitioning.o \
		./build_ins/Read_HDF5Mesh_HOPR.o \
		./build_ins/Read_SpecMesh.o \
		./build_ins/Read_GMSH.o \
		./build_ins/ReadMeshFile.o


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
		./build_ch/Read_GMSH.o \
		./build_ch/ReadMeshFile.o

LIBOBJS_MU =  ./build_mu/ConnectivityClass.o \
		./build_mu/MeshTypes.o \
		./build_mu/WallDistance.o \
		./build_mu/StorageClass.o \
		./build_mu/FacePatchClass.o \
		./build_mu/TransfiniteMaps3D.o \
		./build_mu/MappedGeometry.o \
		./build_mu/NodeClass.o \
		./build_mu/HexElementConnectivityDefinitions.o \
		./build_mu/FaceClass.o \
		./build_mu/HexElementClass.o \
		./build_mu/ZoneClass.o \
		./build_mu/HexMesh.o \
		./build_mu/MeshPartitioning.o \
		./build_mu/Read_HDF5Mesh_HOPR.o \
		./build_mu/Read_SpecMesh.o \
		./build_mu/Read_GMSH.o \
		./build_mu/ReadMeshFile.o

LIB = mesh

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:=all

all: lib$(LIB)_ns.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_mu.a

./build_ns/METISPartitioning.o: ./METISPartitioning.f90		# METIS calls require this, and only this flags... so it is separated
	@echo
	@echo $<
	$(FC) -cpp -O2 -ffree-line-length-0 -DNAVIERSTOKES -I./include_ns $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ns -c $< -o $@ $(MOD_FLAG) ./include_ns

./build_ins/METISPartitioning.o: ./METISPartitioning.f90		# METIS calls require this, and only this flags... so it is separated
	@echo
	@echo $<
	$(FC) -cpp -O2 -ffree-line-length-0 -DINCNS -I./include_ins $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ins -c $< -o $@ $(MOD_FLAG) ./include_ins

./build_ch/METISPartitioning.o: ./METISPartitioning.f90		# METIS calls require this, and only this flags... so it is separated
	@echo
	@echo $<
	$(FC) -cpp -O2 -ffree-line-length-0 -DCAHNHILLIARD -I./include_ch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_ch -c $< -o $@ $(MOD_FLAG) ./include_ch

./build_mu/METISPartitioning.o: ./METISPartitioning.f90		# METIS calls require this, and only this flags... so it is separated
	@echo
	@echo $<
	$(FC) -cpp -O2 -ffree-line-length-0 -DNAVIERSTOKES -DCAHNHILLIARD -I./include_mu $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_mu -c $< -o $@ $(MOD_FLAG) ./include_mu

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

./build_mu/%.o: %.f90
	@echo
	@echo $<
	$(FC) $(FFLAGS) -DMULTIPHASE -DCAHNHILLIARD -I./include_mu $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_mu -c $< -o $@ $(MOD_FLAG) ./include_mu

lib$(LIB)_ns.a: header mkdirs $(LIBOBJS_NS) ./build_ns/METISPartitioning.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_NS) ./build_ns/METISPartitioning.o
	@$(RANLIB) $@

lib$(LIB)_ins.a: header mkdirs $(LIBOBJS_iNS) ./build_ins/METISPartitioning.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_iNS) ./build_ins/METISPartitioning.o
	@$(RANLIB) $@

lib$(LIB)_ch.a: header mkdirs $(LIBOBJS_CH) ./build_ch/METISPartitioning.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_CH) ./build_ch/METISPartitioning.o
	@$(RANLIB) $@

lib$(LIB)_mu.a: header mkdirs $(LIBOBJS_MU) ./build_mu/METISPartitioning.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_MU) ./build_mu/METISPartitioning.o
	@$(RANLIB) $@

install: all
	cp -p lib$(LIB)_ns.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_ins.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_ch.a $(INSTALL_DIR)/lib/
	cp -p lib$(LIB)_mu.a $(INSTALL_DIR)/lib/
	cp -p ./include_ns/*.mod $(INSTALL_DIR)/include_ns/
	cp -p ./include_ins/*.mod $(INSTALL_DIR)/include_ins/
	cp -p ./include_ch/*.mod $(INSTALL_DIR)/include_ch/
	cp -p ./include_mu/*.mod $(INSTALL_DIR)/include_mu/
	@echo
	@echo

install-mu: lib$(LIB)_mu.a
	cp -p lib$(LIB)_mu.a $(INSTALL_DIR)/lib/
	cp -p ./include_mu/*.mod $(INSTALL_DIR)/include_mu/
	@echo
	@echo

install-ns: lib$(LIB)_ns.a
	cp -p lib$(LIB)_ns.a $(INSTALL_DIR)/lib/
	cp -p ./include_ns/*.mod $(INSTALL_DIR)/include_ns/
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
	@mkdir -p ./include_ns ./include_ins ./include_ch ./include_mu
	@mkdir -p ./build_ns ./build_ins ./build_ch ./build_mu

clean: FORCE
	$(RM) lib$(LIB)_ns.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_mu.a
	$(RM) ./build_ns/*.o ./include_ns/*.mod 
	$(RM) ./build_ins/*.o ./include_ins/*.mod 
	$(RM) ./build_ch/*.o ./include_ch/*.mod
	$(RM) ./build_mu/*.o ./include_mu/*.mod

#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE:
