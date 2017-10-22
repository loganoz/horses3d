#####################################################
# Compiling script for the FTObject library
#####################################################

## Add the source files
file(GLOB_RECURSE ft_src "*.f90")

ADD_LIBRARY(${FTLIB} STATIC ${ft_src})
INSTALL(TARGETS ${FTLIB} ARCHIVE DESTINATION ${LIBDIR})
