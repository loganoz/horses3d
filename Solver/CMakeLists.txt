#
# CMake configuration script for HORSES3D
#
#  By: Andrés Rueda
#

##################################################
# Define the project and the depencies that it has
##################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(horses3d Fortran)

# Version of horses3d
SET(VERSION 0.1)

# Default values for additional features
SET(WITH_MKL FALSE CACHE BOOL "Use Intel Math Kernel library")
SET(WITH_LAPACK FALSE CACHE BOOL "Use LAPACK")
SET(WITH_PETSC FALSE CACHE BOOL "Use PETSC")
SET(WITH_HDF5 FALSE CACHE BOOL "Use HDF5")


# Add our local modlues to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Uncomment if it is required that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

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

############################################################
# Locate additional libraries if specified by the user
############################################################
ADD_DEFINITIONS(-DNAVIERSTOKES)

IF(WITH_MKL)
  # Find MKL
  FIND_PACKAGE(MKL REQUIRED)
  SET(HAS_MKL TRUE)
  SET(HAS_LAPACK TRUE)
  MARK_AS_ADVANCED(HAS_MKL)
  MARK_AS_ADVANCED(HAS_LAPACK)
  
  # Macros:
  ADD_DEFINITIONS(-DHAS_MKL)
  ADD_DEFINITIONS(-DHAS_LAPACK)
  # Set BLAS and LAPACK to point to MKL
  SET(BLAS_LIBRARIES ${MKL_BLAS_LIBRARIES} CACHE FILEPATH "")
  SET(LAPACK_LIBRARIES ${MKL_LAPACK_LIBRARIES} CACHE FILEPATH "")
  IF(MPI_FOUND)
    SET(SCALAPACK_LIBRARIES ${MKL_SCALAPACK_LIBRARIES})
    # Versions >11.1 of MKL include Cluster PARDISO direct solver
    IF(WITH_OpenMP AND MKL_CPARDISO_FOUND)
      SET(HAVE_CPARDISO TRUE)
      MARK_AS_ADVANCED(HAVE_CPARDISO)
    ENDIF()
  ENDIF()
  INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR})
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MKL_Fortran_FLAGS}")
  
  SET(CMAKE_REQUIRED_LIBRARIES 
      ${CMAKE_REQUIRED_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
ENDIF()

IF(WITH_LAPACK AND NOT WITH_MKL)
  FIND_PACKAGE(LAPACK REQUIRED)
  ADD_DEFINITIONS(-DHAS_LAPACK)
  SET(HAS_LAPACK TRUE)
  MARK_AS_ADVANCED(HAS_LAPACK)
  # For some reason FindLAPACK doesn't create the LIBS variables in CACHE
  SET (LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE FILEPATH "LAPACK LIBS" FORCE)
  SET (BLAS_LIBRARIES ${BLAS_LIBRARIES} CACHE FILEPATH "BLAS LIBS" FORCE)
ENDIF()

IF(WITH_PETSC)
  FIND_PACKAGE(PETSC REQUIRED)
  ADD_DEFINITIONS(-DHAS_PETSC)
  SET(HAS_PETSC TRUE)
  MARK_AS_ADVANCED(HAS_PETSC)
ENDIF()

IF(WITH_HDF5)
  FIND_PACKAGE(HDF5 REQUIRED)
  ADD_DEFINITIONS(-DHAS_HDF5)
  SET(HAS_HDF5 TRUE)
  MARK_AS_ADVANCED(HAS_HDF5)
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
  SET (HDF5_C_INCLUDE_DIR ${HDF5_C_INCLUDE_DIR} CACHE FILEPATH "HDF5 Fortran INC" FORCE)
  SET (HDF5_Fortran_COMPILER_EXECUTABLE ${HDF5_Fortran_COMPILER_EXECUTABLE} CACHE FILEPATH "HDF5 Fortran COMPILER" FORCE)
  SET (HDF5_Fortran_LIBRARIES ${HDF5_Fortran_LIBRARIES} CACHE FILEPATH "HDF5 Fortran Libraries" FORCE)
#	 SET(CMAKE_Fortran_COMPILER ${HDF5_Fortran_COMPILER_EXECUTABLE})
ENDIF()

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

# Define the executable(s) name(s)
SET(NSEXE HORSES3D)
SET(HORSES2TEC horses2tecplot)
SET(HORSESGEOM horses.geometry)
SET(HORSESF90 horsesf90)

# Define the library names
SET(FOUNDLIB foundation)
SET(FTLIB ftobject)
SET(IOLIB io)
SET(JACLIB jacobian)
SET(MESHLIB mesh)
SET(MONLIB monitors)
SET(PHYSLIB physics)
SET(SPECLIB spectral)
SET(PROBLIB ProblemFile)
SET(MPIUTILLIB MPIUtilities)
SET(DISCLIB discretization)
SET(TIMEINTLIB timeintegrator)

# Define some directories
SET(SRCDIR ${CMAKE_SOURCE_DIR}/src)
SET(LIBDIR ${CMAKE_SOURCE_DIR}/lib)
SET(BINDIR ${CMAKE_SOURCE_DIR}/bin)
SET(INCDIR ${CMAKE_SOURCE_DIR}/include)

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
ADD_SUBDIRECTORY(${SRCDIR}/libs/physics/navierstokes)
ADD_SUBDIRECTORY(${SRCDIR}/libs/spectral)
ADD_SUBDIRECTORY(${SRCDIR}/libs/problemfile)
ADD_SUBDIRECTORY(${SRCDIR}/libs/mpiutils)
ADD_SUBDIRECTORY(${SRCDIR}/libs/discretization)
ADD_SUBDIRECTORY(${SRCDIR}/libs/timeintegrator)

# The source for the main program and add-ons
ADD_SUBDIRECTORY(${SRCDIR}/NavierStokesSolver)
ADD_SUBDIRECTORY(${SRCDIR}/addons)
ADD_SUBDIRECTORY(${SRCDIR}/addons/horses2tecplot)
ADD_SUBDIRECTORY(${SRCDIR}/addons/extractgeometry)

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)
