
#####################################################
# Compiling script for extractgeometry
#####################################################

## Add the source files
#implicitly (easy way -done here):
file(GLOB_RECURSE extractgeometry_src "*.f90")
#explicitly (hard but recommended way -edit this):
#SET(ConsSolver_src DGSEMclass.f90 Horses...)

ADD_EXECUTABLE(${HORSESGEOM} ${extractgeometry_src})

target_include_directories(${HORSESGEOM} PUBLIC ${INCDIR})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################
TARGET_LINK_LIBRARIES(${HORSESGEOM} ${FOUNDLIB} ${MESHLIB} ${SPECLIB} ${FTLIB} ${MPIUTILLIB})


#####################
# Install the binary
#####################

INSTALL(TARGETS ${HORSESGEOM} RUNTIME DESTINATION ${BINDIR})
