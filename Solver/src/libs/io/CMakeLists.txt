#####################################################
# Compiling script for the IO library
#####################################################

## Add the source files
file(GLOB_RECURSE io_src "*.f90")

ADD_LIBRARY(${IOLIB} STATIC ${io_src})
TARGET_LINK_LIBRARIES(${IOLIB} ${FOUNDLIB} ${FTLIB} ${MPIUTILLIB})

target_include_directories(${IOLIB} PUBLIC ${INCDIR})

INSTALL(TARGETS ${IOLIB} ARCHIVE DESTINATION ${LIBDIR})
