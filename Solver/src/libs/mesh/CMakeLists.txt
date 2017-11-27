
###########################################
# Compiling script for the mesh library
###########################################

# Locally define Fortran compiler...

## Add the source files 
file(GLOB_RECURSE mesh_src "*.f90")

ADD_LIBRARY(${MESHLIB} STATIC ${mesh_src})

TARGET_LINK_LIBRARIES(${MESHLIB} ${FOUNDLIB} ${SPECLIB} ${PHYSLIB} ${FTLIB} 
											${IOLIB} ${MPIUTILLIB})
target_include_directories(${MESHLIB} PUBLIC ${INCDIR})

IF(HAS_HDF5)
	TARGET_LINK_LIBRARIES(${MESHLIB} ${HDF5_Fortran_LIBRARIES} ${ZLIB_LIBRARY_RELEASE})
	target_include_directories(${MESHLIB} PUBLIC ${HDF5_C_INCLUDE_DIR})
ENDIF()

INSTALL(TARGETS ${MESHLIB} ARCHIVE DESTINATION ${LIBDIR})
