# ftobjectlibrary

## Add the source files (one can explicitly or implicitly select the source files)
#implicitly (easy way -done here):
file(GLOB_RECURSE dgsem_src "*.f90")
#explicitly (hard but recommended way -edit this):
#SET(ft_src LegendreAlgorithms.f90 Physics.f90 TimeIntegration.f90
#  DGSEM_class.f90 MainProgram.f90 PolynomialModule.f90 Types.f90
#  Input.f90 Parameters.f90 TheRest.f90 Utilities.f90)

ADD_LIBRARY(${DGSEMLIB} STATIC ${dgsem_src})
TARGET_LINK_LIBRARIES(${DGSEMLIB} ${MESHLIB} ${FOUNDLIB} ${JACLIB})

# Link MKL libraries if specified by the user
IF(HAS_MKL)
	TARGET_LINK_LIBRARIES(${DGSEMLIB} ${MKL_LAPACK_LIBRARIES} ${MKL_BLAS_LIBRARIES} ${MKL_SCALAPACK_LIBRARIES})
ENDIF()

IF(HAS_PETSC)
	TARGET_LINK_LIBRARIES(${DGSEMLIB} ${PETSC_LIBRARIES})
ENDIF()

INSTALL(TARGETS ${DGSEMLIB} ARCHIVE DESTINATION ${LIBDIR})
