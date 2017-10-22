
#####################################################
# Compiling script for the conservative solver
#####################################################

## Add the source files
#implicitly (easy way -done here):
file(GLOB_RECURSE ConsSolver_src "*.f90")
#explicitly (hard but recommended way -edit this):
#SET(ConsSolver_src DGSEMclass.f90 Horses...)

ADD_EXECUTABLE(${NSEXE} ${ConsSolver_src})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# This links NSLITE3D to the ftobjectlibrary
TARGET_LINK_LIBRARIES(${NSEXE} ${FTLIB} ${FOUNDLIB} ${FTLIB} ${IOLIB} ${JACLIB} 
                      ${MESHLIB} ${MONLIB} ${PHYSLIB} ${PLOTLIB} ${SPECLIB} ${PROBLIB})

# Link MKL libraries if specified by the user
IF(HAS_MKL)
	TARGET_LINK_LIBRARIES(${NSEXE} ${MKL_LAPACK_LIBRARIES} ${MKL_BLAS_LIBRARIES} ${MKL_SCALAPACK_LIBRARIES})
ENDIF()

IF(HAS_PETSC)
	TARGET_LINK_LIBRARIES(${NSEXE} ${PETSC_LIBRARIES})
ENDIF()

#####################
# Install the binary
#####################

INSTALL(TARGETS ${NSEXE} RUNTIME DESTINATION ${BINDIR})
