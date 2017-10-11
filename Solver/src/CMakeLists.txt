########################################
# Set up how to compile the source files
########################################

## Add the source files
SET(horses_src HORSES3DMain.f90)

#implicitly file(GLOB_RECURSE dgsem_src "*.f90")

# Define the executable in terms of the source files
ADD_EXECUTABLE(${NSEXE} ${horses_src})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# This links NSLITE3D to the ftobjectlibrary
TARGET_LINK_LIBRARIES(${NSEXE} ${FTLIB} ${FOUNDLIB} ${FTLIB} ${IOLIB} ${JACLIB} 
                      ${MESHLIB} ${MONLIB} ${PHYSLIB} ${PLOTLIB} ${SPECLIB} ${PROBLIB} ${DGSEMLIB})
                      

###################################
# Create horsesf90 compiling script
###################################

CONFIGURE_FILE(horsesf90.in.cmake horsesf90)

#########
# Install
#########

INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/horsesf90 DESTINATION ${BINDIR})

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
