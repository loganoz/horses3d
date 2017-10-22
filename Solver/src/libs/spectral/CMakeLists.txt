#####################################################
# Compiling script for the spectral library
#####################################################

## Add the source files
file(GLOB_RECURSE spec_src "*.f90")

ADD_LIBRARY(${SPECLIB} STATIC ${spec_src})
INSTALL(TARGETS ${SPECLIB} ARCHIVE DESTINATION ${LIBDIR})
