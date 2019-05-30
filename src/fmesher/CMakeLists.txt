option(USE_OPENMP "Use OpenMP" OFF)
option(BUILD_FMESHER_WITH_X "Build fmesher with X11 support" OFF)
option(BUILD_FMESHER_WITH_GSL "Build fmesher with GSL support" OFF)
option(BUILD_FMESHER_FOR_DEBUGGING "Build fmesher with DEBUG activated" OFF)

if (USE_OPENMP)
  FIND_PACKAGE( OpenMP)
  if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif()

if(BUILD_FMESHER_WITH_GSL)
  FIND_PACKAGE(GSL)
  if(GSL_FOUND)
    message("GSL FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GSL_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GSL_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GSL_EXE_LINKER_FLAGS}")
  endif()
else()
  ## Turn off the only use of gsl in fmesher:
  add_definitions(-DNO_SPHERICAL_HARMONICS=1)
endif()

if(BUILD_FMESHER_FOR_DEBUGGING)
  message("DEBUG ACTIVATED")
  add_definitions(-DDEBUG)
  add_definitions(-DFMESHER_DEBUG)
endif()

SET_PROPERTY( SOURCE xtmpl.c x11utils.cc PROPERTY COMPILE_FLAGS -Wno-unused-result )
SET_PROPERTY( SOURCE xtmpl.c x11utils.cc PROPERTY COMPILE_FLAGS -Wno-pedantic )

add_definitions(-DLINUX=1)

if(BUILD_FMESHER_WITH_X)
  set(FMESHER_X_FILES xtmpl.c x11utils.cc)
  set(FMESHER_LIBS ${FMESHER_LIBS} X11)
  add_definitions(-DFMESHER_USE_X)
else()
  ## Turn off using the X11 drawning feature:
  add_definitions(-DFMESHER_NO_X)
endif()

# Define executables #
add_library(fmesher STATIC
  src/predicates.cc
  src/mesh.cc
  src/meshc.cc
  src/vector.cc
  src/ioutils.cc
  src/cmdline.c
  src/basis.cc
  src/locator.cc
  src/trees.cc
  ${FMESHER_X_FILES})


## fmesher program for inla
add_executable(fmesher_inla
  src/fmesher.cc)
add_executable(treetest
  src/treetest.cc)
target_link_libraries(fmesher_inla fmesher ${FMESHER_LIBS} ${GSL_LIBRARIES})
target_link_libraries(treetest fmesher ${FMESHER_LIBS} ${GSL_LIBRARIES})

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc_fmesher
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    COMMAND make -C ${CMAKE_CURRENT_BINARY_DIR}/doc/latex
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

#add_subdirectory(test)
