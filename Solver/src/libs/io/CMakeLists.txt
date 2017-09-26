# ftobjectlibrary

## Add the source files (one can explicitly or implicitly select the source files)
#implicitly (easy way -done here):
file(GLOB_RECURSE io_src "*.f90")
#explicitly (hard but recommended way -edit this):
#SET(ft_src LegendreAlgorithms.f90 Physics.f90 TimeIntegration.f90
#  DGSEM_class.f90 MainProgram.f90 PolynomialModule.f90 Types.f90
#  Input.f90 Parameters.f90 TheRest.f90 Utilities.f90)

ADD_LIBRARY(${IOLIB} STATIC ${io_src})
INSTALL(TARGETS ${IOLIB} ARCHIVE DESTINATION ${LIBDIR})
