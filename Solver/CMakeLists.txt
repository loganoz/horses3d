# CMake project file for nslite3d

##################################################
# Define the project and the depencies that it has
##################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(horses3d Fortran)

# Version of horses3d
SET(VERSION 0.1)

# Add our local modlues to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Uncomment if it is required that Fortran 90 is supported
#IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
#    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
#ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
#OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
#OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)

# This INCLUDE statement executes code that sets the compile flags for DEBUG and
# RELEASE. To modify the flags, edit the file SetFortranFlags.cmake.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)

SET(HORSES_F90FLAGS "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}")
 
# Locate and set parallelization libraries.  There are some CMake peculiarities
# taken care of here, such as the fact that the FindOpenMP routine doesn't know
# about Fortran.
#INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)
# Setup the LAPACK libraries.  This also takes care of peculiarities, such as
# the fact the searching for MKL requires a C compiler, and that the results
# are not stored in the cache. 
#INCLUDE(${CMAKE_MODULE_PATH}/SetUpLAPACK.cmake)

# There is an error in CMAKE with this flag for pgf90.  Unset it
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
IF(FCNAME STREQUAL "pgf90")
    UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")

############################
# Add important definitions
############################
ADD_DEFINITIONS(-D_has_Quad)

############################################################
# Define the actual files and folders that make up the build
############################################################

# Define the executable(s) name(s)... arueda: distribute with SpecMesh3D?
SET(NSEXE HORSES3D)
SET(HORSESF90 horsesf90)

# Define the library names
SET(FOUNDLIB foundation)
SET(FTLIB ftobject)
SET(IOLIB io)
SET(JACLIB jacobian)
SET(MESHLIB mesh)
SET(MONLIB monitors)
SET(PHYSLIB physics)
SET(PLOTLIB plot)
SET(SPECLIB spectral)
SET(PROBLIB ProblemFile)
SET(DGSEMLIB dgsem)

# Define some directories
SET(SRCDIR ${CMAKE_SOURCE_DIR}/src)
SET(LIBDIR ${CMAKE_SOURCE_DIR}/lib)
SET(BINDIR ${CMAKE_SOURCE_DIR}/bin)
SET(INCDIR ${CMAKE_SOURCE_DIR}/include)

#---------------
#To be deprecated:
#SET(SRCFOO ${SRC}/foo)
#SET(SRCBAR ${SRC}/bar)
#---------------

# Have the .mod files placed in the include folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${INCDIR})

# The source for the libraries
ADD_SUBDIRECTORY(${SRCDIR}/libs/foundation)
ADD_SUBDIRECTORY(${SRCDIR}/libs/ftobject)
ADD_SUBDIRECTORY(${SRCDIR}/libs/io)
ADD_SUBDIRECTORY(${SRCDIR}/libs/jacobian)
ADD_SUBDIRECTORY(${SRCDIR}/libs/mesh)
ADD_SUBDIRECTORY(${SRCDIR}/libs/monitors)
ADD_SUBDIRECTORY(${SRCDIR}/libs/physics)
ADD_SUBDIRECTORY(${SRCDIR}/libs/plot)
ADD_SUBDIRECTORY(${SRCDIR}/libs/spectral)
ADD_SUBDIRECTORY(${SRCDIR}/libs/problemfile)
ADD_SUBDIRECTORY(${SRCDIR}/libs/dgsem)

# The source for the NSLITE3D binary
ADD_SUBDIRECTORY(${SRCDIR})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)
