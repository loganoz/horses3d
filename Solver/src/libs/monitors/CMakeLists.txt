
#####################################################
# Compiling script for the monitors library
#####################################################

## Add the source files
file(GLOB_RECURSE mon_src "*.f90")

ADD_LIBRARY(${MONLIB} STATIC ${mon_src})
TARGET_LINK_LIBRARIES(${MONLIB} ${MESHLIB} ${FTLIB})
INSTALL(TARGETS ${MONLIB} ARCHIVE DESTINATION ${LIBDIR})
