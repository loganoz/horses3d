
###########################################
# Compiling script for the mesh library
###########################################

## Add the source files 
file(GLOB_RECURSE mesh_src "*.f90")

ADD_LIBRARY(${MESHLIB} STATIC ${mesh_src})

TARGET_LINK_LIBRARIES(${MESHLIB} ${SPECLIB} ${PHYSLIB} ${FTLIB} ${IOLIB})

INSTALL(TARGETS ${MESHLIB} ARCHIVE DESTINATION ${LIBDIR})
