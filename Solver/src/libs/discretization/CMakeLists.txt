
#################################################
# Compiling script for the discretization library
#################################################

## Add the source files 
file(GLOB_RECURSE disc_src "*.f90")

ADD_LIBRARY(${DISCLIB} STATIC ${disc_src})

TARGET_LINK_LIBRARIES(${DISCLIB} ${FOUNDLIB} ${PHYSLIB} ${MESHLIB} ${MPIUTILLIB} ${FTLIB} ${IOLIB})
target_include_directories(${DISCLIB} PUBLIC ${INCDIR})

INSTALL(TARGETS ${DISCLIB} ARCHIVE DESTINATION ${LIBDIR})
