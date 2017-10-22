#####################################################
# Compiling script for the IO library
#####################################################

## Add the source files
file(GLOB_RECURSE io_src "*.f90")

ADD_LIBRARY(${IOLIB} STATIC ${io_src})
INSTALL(TARGETS ${IOLIB} ARCHIVE DESTINATION ${LIBDIR})
