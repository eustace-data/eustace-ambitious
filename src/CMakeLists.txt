cmake_minimum_required (VERSION 2.8)

project(ambitious)

option(USE_OPENMP "Use OpenMP" OFF)
option(ON_CEMS "Running on CEMS" OFF)
option(ACTIVATE_DEBUG "Build DEBUG activated" OFF)
option(ACTIVATE_DEBUGINFO "Build with -g -O0" OFF)
set( EIGEN_PATH "" CACHE PATH "Path to Eigen3 (set empty to use system)")
set( BOOST_PATH "" CACHE PATH "Boost library path (leave empty for system version" )
set( Boost_NO_SYSTEM_PATHS OFF CACHE BOOL "Do not search system for Boost" )

if (ON_CEMS)
  include_directories("/apps/libs/netCDF/intel17/4.6.1/include")
  link_directories("/apps/libs/netCDF/intel17/4.6.1/lib/")

  include_directories("/apps/libs/PHDF5/intel17/1.10.0_patch1/include")
  link_directories("/apps/libs/PHDF5/intel17/1.10.0_patch1/lib/")

  set(PNETCDF_INCLUDE_PATH "/apps/libs/pnetcdf/intel17/1.10.0/include/")
  link_directories("/apps/libs/pnetcdf/intel17/1.10.0/lib/")

  if(NOT CMAKE_COMPILER_IS_GNUCXX)
    ## This test only works because the cems system uses an old Intel compiler
    add_definitions(-DOLDCXX=1)
    message("Compiler without full c++11 support detected. Defining OLDCXX.")
  endif(NOT CMAKE_COMPILER_IS_GNUCXX)
endif()

# Windows configuration parameters
if (MSVC)
  option(EUSTACE_MSVC_RUNTIME_STATIC "Use static (non-DLL) version of C runtime libraries" ON)
endif(MSVC)

# Optionally use static libraries on windows
if (EUSTACE_MSVC_RUNTIME_STATIC)
  foreach(flag_var
        CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
        CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO 
	CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
        CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO)
    if(${flag_var} MATCHES "/MD")
      string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
    endif(${flag_var} MATCHES "/MD")
  endforeach(flag_var)
endif(EUSTACE_MSVC_RUNTIME_STATIC)


FIND_PACKAGE( OpenMP)
if(USE_OPENMP)
  if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif()

FIND_PACKAGE( MPI)
if(MPI_FOUND)
  message("MPI FOUND")
  link_directories("${PROJECT_SOURCE_DIR}/../linker/")
endif()

if(ACTIVATE_DEBUG)
  message("DEBUG ACTIVATED")
  add_definitions(-DDEBUG)
endif()

if(ACTIVATE_DEBUGINFO)
  message("DEBUGINFO ACTIVATED")
  set(DEBUGINFO_FLAGS "-g")
  set(OPTIMISATION_FLAGS "-O0")
else()
  set(DEBUGINFO_FLAGS "")
  set(OPTIMISATION_FLAGS "-O2")
endif()


########## Find boost #################
if(NOT BOOST_PATH STREQUAL "")
  message("Setting BOOST_ROOT")
  set(BOOST_ROOT ${BOOST_PATH})
endif()
FIND_PACKAGE(Boost COMPONENTS program_options system REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
########## Find boost #################


# Setup Eigen #
set( EIGEN_TARGET "Eigen3::Eigen")
if(NOT EIGEN_PATH STREQUAL "")
  message("Adding Eigen to cmake include path")
  list(APPEND CMAKE_INCLUDE_PATH ${EIGEN_PATH})
endif()
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
###############

# Setup paths for the eustace and fmesher modules
link_directories("${PROJECT_BINARY_DIR}/eustace")
include_directories("${PROJECT_SOURCE_DIR}/fmesher/src/")
link_directories("${PROJECT_BINARY_DIR}/fmesher")

# Setup common C++ compiler options
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${DEBUGINFO_FLAGS} ${OPTIMISATION_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DEBUGINFO_FLAGS} ${OPTIMISATION_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpmath=sse -msse2")


include_directories ("${PROJECT_SOURCE_DIR}")
add_subdirectory(eustace)
add_subdirectory(fmesher)
add_subdirectory(ambitious)
