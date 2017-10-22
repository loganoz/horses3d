#####################################################
# Compiling script for the foundation library
#####################################################

## Add the source files
file(GLOB_RECURSE found_src "*.f90")

ADD_LIBRARY(${FOUNDLIB} STATIC ${found_src})
TARGET_LINK_LIBRARIES(${FOUNDLIB} ${FTLIB})

# Copy Includes.h to include directory
file(COPY Includes.h
        DESTINATION ${INCDIR})

# Link LAPACK libraries if specified by the user
IF(HAS_LAPACK)
  TARGET_LINK_LIBRARIES(${FOUNDLIB} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
ENDIF()

INSTALL(TARGETS ${FOUNDLIB} ARCHIVE DESTINATION ${LIBDIR})
