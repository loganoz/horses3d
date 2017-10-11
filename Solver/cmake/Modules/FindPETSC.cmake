# Script for finding the PETSC library
#	  (very simplified version based on https://github.com/jedbrown/cmake-modules/blob/master/FindPETSc.cmake)
# By Andrés Rueda
#
#
###############
# Get PETSc DIR
###############

find_path (PETSC_DIR include/petsc.h
  HINTS ENV PETSC_DIR
  PATHS
  ${CMAKE_SYSTEM_PREFIX_PATH}
  # Debian paths
  /usr/lib/petscdir/3.5.1 /usr/lib/petscdir/3.5
  /usr/lib/petscdir/3.4.2 /usr/lib/petscdir/3.4
  /usr/lib/petscdir/3.3 /usr/lib/petscdir/3.2 /usr/lib/petscdir/3.1
  /usr/lib/petscdir/3.0.0 /usr/lib/petscdir/2.3.3 /usr/lib/petscdir/2.3.2
  # MacPorts path
  /opt/local/lib/petsc
  DOC "PETSc Directory")
  
###################################
# Get architecture (directory name)
# ... looks for env variable... If there's not one defined, I give some hints
###################################

# If not defined by the user, assign "NOT-FOUND"
if (NOT PETSC_ARCH)
  set (PETSC_ARCH "NOT-FOUND" CACHE STRING "PETSc build architecture")
endif (NOT PETSC_ARCH)

# If "NOT-FOUND", then try to find it... 
if (PETSC_DIR AND PETSC_ARCH STREQUAL "NOT-FOUND")
  set (_petsc_arches
    $ENV{PETSC_ARCH}                   # If set, use environment variable first
    linux-gnu-c-debug linux-gnu-c-opt  # Debian defaults
    linux-gnu-f-debug linux-gnu-f-opt  # Debian defaults
    linux-gnu-mpich linux-ifort-mpich  # My paths...
    x86_64-unknown-linux-gnu i386-unknown-linux-gnu)
    
  set (petscconf "NOTFOUND" CACHE FILEPATH "Cleared" FORCE)
  foreach (arch ${_petsc_arches})
    if (PETSC_ARCH STREQUAL "NOT-FOUND")
      find_path (petscconf petscconf.h
        HINTS ${PETSC_DIR}
        PATH_SUFFIXES ${arch}/include bmake/${arch}
        NO_DEFAULT_PATH)
      if (petscconf)
        set (PETSC_ARCH "${arch}" CACHE STRING "PETSc build architecture" FORCE)
      endif (petscconf)
    endif (PETSC_ARCH STREQUAL "NOT-FOUND")
  endforeach (arch)
  set (petscconf "NOTFOUND" CACHE INTERNAL "Scratch variable" FORCE)
endif (PETSC_DIR AND PETSC_ARCH STREQUAL "NOT-FOUND")

# Now if that didn't find it, then the user will realize they have to specify it... unless the installation is in the petsc_dir


###################################
# Determine whether the PETSc layout is old-style (through 2.3.3) or
# new-style (3.0.0)
###################################

if (EXISTS "${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/petscvariables") # > 3.5
  
#  set (PETSC_INCLUDES "${PETSC_DIR}/${PETSC_ARCH}/include ${PETSC_DIR}/include" CACHE STRING "PETSc include path" FORCE)
  set (PETSC_LIBRARIES "${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.so" CACHE FILEPATH "PETSc libraries" FORCE)
  
  INCLUDE_DIRECTORIES(${PETSC_DIR}/${PETSC_ARCH}/include)
  INCLUDE_DIRECTORIES(${PETSC_DIR}/include)
  
  set (PETSC_FOUND Yes)
  
elseif (EXISTS "${PETSC_DIR}/bmake/${PETSC_ARCH}/petscconf.h") # <= 2.3.3
#  set (PETSC_INCLUDES "${PETSC_DIR}/bmake/${PETSC_ARCH}/include ${PETSC_DIR}/include" CACHE STRING "PETSc include path" FORCE)
  set (PETSC_LIBRARIES "${PETSC_DIR}/bmake/${PETSC_ARCH}/lib/libpetsc.so" CACHE FILEPATH "PETSc libraries" FORCE)
  
  INCLUDE_DIRECTORIES(${PETSC_DIR}/bmake/${PETSC_ARCH}/include)
  INCLUDE_DIRECTORIES(${PETSC_DIR}/include)
  
  set (PETSC_FOUND Yes)
elseif (PETSC_DIR)
  message (SEND_ERROR "The pair PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} do not specify a valid PETSc installation")
endif ()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSC DEFAULT_MSG PETSC_LIBRARIES) #PETSC_INCLUDES
