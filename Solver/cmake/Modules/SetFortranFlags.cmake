######################################################
# Determine and set the Fortran compiler flags we want 
######################################################

####################################################################
# Make sure that the default build type is RELEASE if not specified.
####################################################################
INCLUDE(${CMAKE_MODULE_PATH}/SetCompileFlag.cmake)

# Make sure the build type is uppercase
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

IF(BT STREQUAL "RELEASE")
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
ELSEIF(BT STREQUAL "DEBUG")
    SET (CMAKE_BUILD_TYPE DEBUG CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
#ELSEIF(BT STREQUAL "TESTING") #arueda: outcommented
#    SET (CMAKE_BUILD_TYPE TESTING CACHE STRING
#      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
#      FORCE)
ELSEIF(NOT BT)
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
    MESSAGE(STATUS "CMAKE_BUILD_TYPE not given, defaulting to RELEASE")
ELSE()
    MESSAGE(FATAL_ERROR "CMAKE_BUILD_TYPE not valid, choices are DEBUG, or RELEASE") #, or TESTING
ENDIF(BT STREQUAL "RELEASE")

#########################################################
# If the compiler flags have already been set, return now
#########################################################

IF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_TESTING AND CMAKE_Fortran_FLAGS_DEBUG)
    RETURN ()
ENDIF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_TESTING AND CMAKE_Fortran_FLAGS_DEBUG)

########################################################################
# Determine the appropriate flags for this compiler for each build type.
# For each option type, a list of possible flags is given that work
# for various compilers.  The first flag that works is chosen.
# If none of the flags work, nothing is added (unless the REQUIRED 
# flag is given in the call).  This way unknown compiles are supported.
#######################################################################

#####################
### GENERAL FLAGS ###
#####################

## Don't add underscores in symbols for C-compatability ... arueda: Will this be needed?
#SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
#                 Fortran "-fno-underscoring")


## Following flags Optimize for the host's architecture
##   arueda: Very good, but then it's highly likely that code on dropbox won't work
# There is some bug where -march=native doesn't work on Mac
IF(APPLE)
    SET(GNUNATIVE "-mtune=native")
ELSE()
    SET(GNUNATIVE "-march=native")
ENDIF()

# Optimize for the host's architecture
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-xHost"        # Intel
                         "/QxHost"       # Intel Windows
                         ${GNUNATIVE}    # GNU
                         "-ta=host"      # Portland Group
                )

## C-preprocessor flag
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran REQUIRED "-cpp"              # Intel/GNU
                                  "/cpp"              # Intel Windows
                                  "-qsuffix=cpp=f90"  # IBM...
                                  #" "         # Portland Group?
                )

# Change the interpretation of backslashes in string literals from a single backslash character to “C-style” escape characters ## DEBUG flag??
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-fbackslash"  # GNU
                         "-assume bscc" # Intel
                         "-qescape"     # IBM...
                )

## fPIC    # This is needed in order to compile the ProblemFile.so as dynamic library (change??)
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran REQUIRED "-fPIC"       # Intel/GNU
                )

## Free line length to avoid "line truncated" error
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-ffree-line-length-0"       # GNU
                )

###################
### DEBUG FLAGS ###
###################

# NOTE: debugging symbols (-g or /debug:full) are already on by default

# Disable optimizations
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran REQUIRED "-O0" # All compilers not on Windows
                                  "/Od" # Intel Windows
                )

# Turn on all warnings 
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-warn all" # Intel
                         "/warn:all" # Intel Windows
                         "-Wall"     # GNU
                                     # Portland Group (on by default)
                )
# But please don't tell me about the unused variables
#SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
#                 Fortran "-Wno-unused-variable"     # GNU
#                         "-warn all [no]unused"    # Intel
#                )

# Traceback
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-traceback"   # Intel/Portland Group
                         "/traceback"   # Intel Windows
                         "-fbacktrace"  # GNU (gfortran)
                         "-ftrace=full" # GNU (g95)
                )
# Check all
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran REQUIRED "-check all"  # Intel
                                  "/check:all"  # Intel Windows
                                  "-fcheck=all" # GNU (New style)
                                  "-C=all"      # nag95
                                  "xcheck=%all" # Sun F95
                                  "--chk"       # If95
                )

# Warn about some not initiallized variables
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-finit-real=snan" # GNU
                          # Intel option?
                )

#
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-ffpe-trap=invalid,zero,overflow" # GNU
                          # Intel option?
                )

#####################
### TESTING FLAGS ### #arueda: outcommented
#####################

# Optimizations
#SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_TESTING "${CMAKE_Fortran_FLAGS_TESTING}"
#                 Fortran REQUIRED "-O2" # All compilers not on Windows
#                                  "/O2" # Intel Windows
#                )

#####################
### RELEASE FLAGS ###
#####################

# NOTE: agressive optimizations (-O3) are already turned on by default

#    -fopenmp

# Unroll loops # arueda: useful?
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-funroll-loops" # GNU
                         "-unroll"        # Intel
                         "/unroll"        # Intel Windows
                         "-Munroll"       # Portland Group
                )

# Inline functions
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-inline"            # Intel
                         "/Qinline"           # Intel Windows
                         "-finline-functions" # GNU
                         "-Minline"           # Portland Group
                )

# Interprocedural (link-time) optimizations  #arueda: outcommented because was generating problems
#SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#                 Fortran "-ipo"     # Intel
#                         "/Qipo"    # Intel Windows
#                         "-flto"    # GNU
#                         "-Mipa"    # Portland Group
#                )



# Single-file optimizations
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-ip"  # Intel
                         "/Qip" # Intel Windows
                )

# Vectorize code
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-qopt-report0"                                   # Intel
                         "/qopt-report0"                                  # Intel Windows
                         "-Mvect"                                         # Portland Group
                         "-ftree-vectorize  -ftree-vectorizer-verbose=0"  # GNU ... Already included in -O3
                )
# OpenMP flags:
IF(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
    FIND_PACKAGE(OpenMP REQUIRED)
    # Here I don't use SET_COMPILE_FLAG function since some bug where OpenMP_Fortran_FLAGS is used (anyway cmake looks for the appropriate flag)
    SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_Fortran_FLAGS}")
    IF (NOT OpenMP_Fortran_FLAGS)
        MESSAGE (FATAL_ERROR "Fortran compiler does not support OpenMP")
    ENDIF (NOT OpenMP_Fortran_FLAGS)
    
    # For IBM -qsmp=omp
    
ENDIF(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
