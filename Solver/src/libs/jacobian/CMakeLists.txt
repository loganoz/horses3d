#####################################################
# Compiling script for the Jacobian library
#####################################################

## Add the source files 
file(GLOB_RECURSE jac_src "*.f90")

ADD_LIBRARY(${JACLIB} STATIC ${jac_src})
TARGET_LINK_LIBRARIES(${JACLIB} ${MESHLIB})

# Link MKL libraries if specified by the user
IF(HAS_MKL)
	TARGET_LINK_LIBRARIES(${JACLIB} ${MKL_LAPACK_LIBRARIES} ${MKL_BLAS_LIBRARIES} ${MKL_SCALAPACK_LIBRARIES})
ENDIF()

INSTALL(TARGETS ${JACLIB} ARCHIVE DESTINATION ${LIBDIR})
