
#####################################################
# Compiling script for the monitors library
#####################################################

## Add the source files
file(GLOB_RECURSE mon_src "*.f90")

ADD_LIBRARY(${MONLIB} STATIC ${mon_src})

target_include_directories(${MONLIB} PUBLIC ${INCDIR})
TARGET_LINK_LIBRARIES(${MONLIB} ${FOUNDLIB} ${FTLIB} ${MESHLIB} ${IOLIB})

INSTALL(TARGETS ${MONLIB} ARCHIVE DESTINATION ${LIBDIR})
