#####################################################
# Compiling script for the Physics library
#####################################################

## Add the source files
SET(phys_src BoundaryConditions.f90 FluidData.f90 ManufacturedSolutions.f90 Physics.f90)

ADD_LIBRARY(${PHYSLIB} STATIC ${phys_src})

TARGET_LINK_LIBRARIES(${PHYSLIB} ${FOUNDLIB} ${FTLIB} ${MESHLIB} ${IOLIB})
target_include_directories(${PHYSLIB} PUBLIC ${INCDIR})

INSTALL(TARGETS ${PHYSLIB} ARCHIVE DESTINATION ${LIBDIR})
