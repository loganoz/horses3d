#####################################################
# Compiling script for the spectral library
#####################################################

## Add the source files
file(GLOB_RECURSE spec_src "*.f90")

ADD_LIBRARY(${SPECLIB} STATIC ${spec_src})

target_include_directories(${SPECLIB} PUBLIC ${INCDIR})

TARGET_LINK_LIBRARIES(${SPECLIB} ${FOUNDLIB})

INSTALL(TARGETS ${SPECLIB} ARCHIVE DESTINATION ${LIBDIR})
