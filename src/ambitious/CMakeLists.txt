# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(ambitious_doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    COMMAND make -C ${CMAKE_CURRENT_BINARY_DIR}/doc/latex
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
endif(DOXYGEN_FOUND)

option(USE_PNETCDF "Use pnetcdf" ON)

set(MYNETCDF_LIBRARIES "netcdf")
if (MPI_FOUND)
  #  set(MYNETCDF_LIBRARIES ${MYNETCDF_LIBRARIES} pnetcdf "${MPI_CXX_LIBRARIES}")
  
  #    add_definitions(-DOMPI_SKIP_MPICXX)
  
  include_directories(${MPI_CXX_INCLUDE_PATH})

  message("MPI_CXX_INCLUDE_PATH: ${MPI_CXX_INCLUDE_PATH}")
  message("MPI_CXX_LIBRARIES: ${MPI_CXX_LIBRARIES}")
  
  if(USE_PNETCDF)
    message("PNETCDF ACTIVATED")
    add_definitions(-DUSE_PNETCDF)
    set(MYNETCDF_LIBRARIES ${MYNETCDF_LIBRARIES} pnetcdf "${MPI_CXX_LIBRARIES}")

    include_directories(${PNETCDF_INCLUDE_PATH})
#    link_libraries(${PNETCDF_CXX_LIBRARIES})

    message("PNETCDF_INCLUDE_PATH: ${PNETCDF_INCLUDE_PATH}")
    message("PNETCDF_CXX_LIBRARIES: ${PNETCDF_CXX_LIBRARIES}")

  endif(USE_PNETCDF)
endif(MPI_FOUND)

add_library(ambitious
  blocks/blocks.cpp
  timer/timer.cpp
  po_tool/po_tool.cpp
  po/po.cpp
  meshes/meshes.cpp
  ncdf/ncdf_helper.cpp
  observe/observe.cpp
  versions/versions.cpp
  debuglog/debuglog.cpp)
target_link_libraries(ambitious PUBLIC ${EIGEN_TARGET})

add_subdirectory(informal)
add_subdirectory(test)
