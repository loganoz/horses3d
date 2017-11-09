#####################################################
# Compiling script for the plot library
#####################################################

## Add the source files
file(GLOB_RECURSE plot_src "*.f90")

ADD_LIBRARY(${PLOTLIB} STATIC ${plot_src})

TARGET_LINK_LIBRARIES(${PLOTLIB} ${FOUNDLIB} ${PHYSLIB} ${SPECLIB}) 
INSTALL(TARGETS ${PLOTLIB} ARCHIVE DESTINATION ${LIBDIR})
