
################################################
# Compiling script for the mpi utilities library
################################################

# Locally define Fortran compiler...

## Add the source files 
file(GLOB_RECURSE mpi_src "*.f90")

ADD_LIBRARY(${MPIUTILLIB} STATIC ${mpi_src})

TARGET_LINK_LIBRARIES(${MPIUTILLIB} ${FOUNDLIB})
target_include_directories(${MPIUTILLIB} PUBLIC ${INCDIR})

INSTALL(TARGETS ${MPIUTILLIB} ARCHIVE DESTINATION ${LIBDIR})
