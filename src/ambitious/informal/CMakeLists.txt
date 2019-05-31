add_executable(cal
  cal.cpp)
add_executable(modelscan
  modelscan.cpp)
add_executable(build_bilevel_meshes
  build_bilevel_meshes.cpp)
add_executable(qtool_timing
  qtool_timing.cpp)
add_executable(commandline_demo
  commandline_demo.cpp)
add_executable(read_one_day_test
  read_one_day_test.cpp)
add_executable(meta_data_in_input
  meta_data_in_input.cpp)
add_executable(obs
  obs.cpp)

set(COMMON_LIBRARIES_SUBSET ambitious ${Boost_LIBRARIES} fmesher eustace)
set(COMMON_LIBRARIES ambitious ${Boost_LIBRARIES} fmesher eustace ${MYNETCDF_LIBRARIES} ${MPI_CXX_LIBRARIES})
message(${MYNETCDF_LIBRARIES})
target_link_libraries(cal ${COMMON_LIBRARIES})
target_link_libraries(modelscan ${COMMON_LIBRARIES})
target_link_libraries(build_bilevel_meshes ${COMMON_LIBRARIES})
target_link_libraries(qtool_timing ${COMMON_LIBRARIES})
target_link_libraries(commandline_demo ${COMMON_LIBRARIES})
target_link_libraries(read_one_day_test ${COMMON_LIBRARIES})
target_link_libraries(meta_data_in_input ${COMMON_LIBRARIES})
target_link_libraries(obs ${COMMON_LIBRARIES})
