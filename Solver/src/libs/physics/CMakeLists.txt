#####################################################
# Compiling script for the Physics library
#####################################################

## Add the source files
file(GLOB_RECURSE phys_src "*.f90")

ADD_LIBRARY(${PHYSLIB} STATIC ${phys_src})

TARGET_LINK_LIBRARIES(${PHYSLIB} ${FOUNDLIB} ${FTLIB} ${MESHLIB} ${IOLIB} ${MPIUTILLIB})
target_include_directories(${PHYSLIB} PUBLIC ${INCDIR})

INSTALL(TARGETS ${PHYSLIB} ARCHIVE DESTINATION ${LIBDIR})
