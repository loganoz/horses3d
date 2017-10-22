
#####################################################
# Compiling script for the conservative solver
#####################################################

## Add the source files
#implicitly (easy way -done here):
file(GLOB_RECURSE horses2tec_src "*.f90")
#explicitly (hard but recommended way -edit this):
#SET(ConsSolver_src DGSEMclass.f90 Horses...)

ADD_EXECUTABLE(${HORSES2TEC} ${horses2tec_src})

target_include_directories(${HORSES2TEC} PUBLIC ${INCDIR})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################
TARGET_LINK_LIBRARIES(${HORSES2TEC} ${FOUNDLIB} ${MESHLIB} ${IOLIB})


#####################
# Install the binary
#####################

INSTALL(TARGETS ${HORSES2TEC} RUNTIME DESTINATION ${BINDIR})
