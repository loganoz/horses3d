########################################
# Set up how to compile the source files
########################################

## Add the source files (one can explicitly or implicitly select the source files)
#implicitly (easy way -done here):
#file(GLOB_RECURSE dgsem_src "*.f90")
#explicitly (hard but recommended way -edit this):
SET(horses_src HORSES3DMain.f90)
# Physics.f90 TimeIntegration.f90
#  DGSEM_class.f90 MainProgram.f90 PolynomialModule.f90 Types.f90
#  Input.f90 Parameters.f90 TheRest.f90 Utilities.f90)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${NSEXE} ${horses_src})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# This links NSLITE3D to the ftobjectlibrary
TARGET_LINK_LIBRARIES(${NSEXE} ${FTLIB} ${FOUNDLIB} ${FTLIB} ${IOLIB} ${JACLIB} 
                      ${MESHLIB} ${MONLIB} ${PHYSLIB} ${PLOTLIB} ${SPECLIB} ${PROBLIB} ${DGSEMLIB})

CONFIGURE_FILE(horsesf90.in.cmake horsesf90)
  						 
INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/horsesf90 DESTINATION ${BINDIR})
# Uncomment if you need to link to BLAS and LAPACK - arueda: not yet
#TARGET_LINK_LIBRARIES(${FOOEXE} ${BLAS_LIBRARIES}
#                                ${LAPACK_LIBRARIES}
#                                ${CMAKE_THREAD_LIBS_INIT})

# Uncomment if you have parallization
#IF(USE_OPENMP)
#    SET_TARGET_PROPERTIES(${FOOEXE} PROPERTIES
#                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
#                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
#ELSEIF(USE_MPI)
#    SET_TARGET_PROPERTIES(${FOOEXE} PROPERTIES
#                          COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
#                          LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
#    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
#    TARGET_LINK_LIBRARIES(${FOOEXE} ${MPI_Fortran_LIBRARIES})
#ENDIF(USE_OPENMP)

#####################################
# Tell how to install this executable
#####################################

#IF(WIN32)
#    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
#ELSE()
#    SET(CMAKE_INSTALL_PREFIX /usr/local)
#ENDIF(WIN32)
INSTALL(TARGETS ${NSEXE} RUNTIME DESTINATION ${BINDIR})
