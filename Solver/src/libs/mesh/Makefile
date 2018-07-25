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

LIBOBJS_iNSCH = ./build_insch/ConnectivityClass.o \
		./build_insch/MeshTypes.o \
		./build_insch/WallDistance.o \
		./build_insch/StorageClass.o \
		./build_insch/FacePatchClass.o \
		./build_insch/TransfiniteMaps3D.o \
		./build_insch/MappedGeometry.o \
		./build_insch/NodeClass.o \
		./build_insch/HexElementConnectivityDefinitions.o \
		./build_insch/FaceClass.o \
		./build_insch/HexElementClass.o \
		./build_insch/ZoneClass.o \
		./build_insch/HexMesh.o \
		./build_insch/MeshPartitioning.o \
		./build_insch/Read_HDF5Mesh_HOPR.o \
		./build_insch/Read_SpecMesh.o \
		./build_insch/ReadMeshFile.o



LIB = mesh

############################################
#	MAIN COMPILATION RULES 		   #
############################################
.DEFAULT_GOAL:=all

all: lib$(LIB)_ns.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_nsch.a lib$(LIB)_insch.a

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

./build_nsch/METISPartitioning.o: ./METISPartitioning.f90		# METIS calls require this, and only this flags... so it is separated
	@echo
	@echo $<
	$(FC) -cpp -O2 -ffree-line-length-0 -DNAVIERSTOKES -DCAHNHILLIARD -I./include_nsch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_nsch -c $< -o $@ $(MOD_FLAG) ./include_nsch

./build_insch/METISPartitioning.o: ./METISPartitioning.f90		# METIS calls require this, and only this flags... so it is separated
	@echo
	@echo $<
	$(FC) -cpp -O2 -ffree-line-length-0 -DINSCH -DCAHNHILLIARD -I./include_insch $(MACROS) $(EXTLIB_INC) $(INCLUDE) $(INCLUDE)_insch -c $< -o $@ $(MOD_FLAG) ./include_insch

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

lib$(LIB)_nsch.a: header mkdirs $(LIBOBJS_NSCH) ./build_nsch/METISPartitioning.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_NSCH) ./build_nsch/METISPartitioning.o
	@$(RANLIB) $@

lib$(LIB)_insch.a: header mkdirs $(LIBOBJS_iNSCH) ./build_insch/METISPartitioning.o
	@echo
	@echo "---------------------------------"
	@echo "| Linking library" $@ "|"
	@echo "---------------------------------"
	@$(RM) $@
	$(AR) $(ARFLAGS) $@ $(LIBOBJS_iNSCH) ./build_insch/METISPartitioning.o
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
	@echo ">>   Building Mesh        "
	@echo "================================"
	@echo
	@echo
mkdirs: FORCE
	@mkdir -p ./include_ns ./include_ins ./include_ch ./include_nsch ./include_insch
	@mkdir -p ./build_ns ./build_ins ./build_ch ./build_nsch ./build_insch

clean: FORCE
	$(RM) lib$(LIB)_ns.a lib$(LIB)_ins.a lib$(LIB)_ch.a lib$(LIB)_nsch.a
	$(RM) ./build_ns/*.o ./include_ns/*.mod 
	$(RM) ./build_ins/*.o ./include_ins/*.mod 
	$(RM) ./build_ch/*.o ./include_ch/*.mod
	$(RM) ./build_nsch/*.o ./include_nsch/*.mod
	$(RM) ./build_insch/*.o ./include_insch/*.mod

#############################################
#	Dummy procedure to force a rule     #
#############################################
FORCE:
