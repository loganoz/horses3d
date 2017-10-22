
###########################################
# Compiling script for Problemfile library
###########################################

file(GLOB_RECURSE dyn_src "*.f90")

ADD_LIBRARY(${PROBLIB} SHARED ${dyn_src})

TARGET_LINK_LIBRARIES(${PROBLIB} ${FTLIB} ${FOUNDLIB} ${MESHLIB})

SET_TARGET_PROPERTIES(${PROBLIB} PROPERTIES PREFIX "")
SET_TARGET_PROPERTIES(${PROBLIB} PROPERTIES LINKER_LANGUAGE Fortran)

INSTALL(TARGETS ${PROBLIB} LIBRARY DESTINATION ${LIBDIR}
                           RUNTIME DESTINATION ${LIBDIR})
