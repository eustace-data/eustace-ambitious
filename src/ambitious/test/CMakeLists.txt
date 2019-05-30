include_directories(/group_workspaces/cems2/eustace/software/cppunit/include)
link_directories(/group_workspaces/cems2/eustace/software/cppunit/lib)

## add_definitions(-DEUSTACE_TEST_DATA_DIRECTORY="${CMAKE_CURRENT_SOURCE_DIR}/data")
#### add_definitions(-fprofile-arcs -ftest-coverage)

add_executable(test_NCDFScaling
  main.cpp
  TestNCDFScaling.cpp)
add_executable(test_blocks
  main.cpp
  TestBlocks.cpp)
add_executable(test_poq
  main.cpp
  TestPoq.cpp)
add_executable(test_timer
  main.cpp
  TestTimer.cpp)
add_executable(test_bidirmap
  main.cpp
  TestBidirmap.cpp)
add_executable(test_qtool
  main.cpp
  TestQtool.cpp)
add_executable(test_po_tool
  main.cpp
  TestPoTool.cpp)
add_executable(test_meshes
  main.cpp
  TestMeshes.cpp)

add_executable(test_all
  main.cpp
  TestNCDFScaling.cpp
  TestBlocks.cpp
  TestPoq.cpp
  TestTimer.cpp
  TestBidirmap.cpp
  TestMeshes.cpp
  TestQtool.cpp
  TestPoTool.cpp)

set(MYLIBRARIES cppunit ambitious ${Boost_LIBRARIES} fmesher eustace ${MYNETCDF_LIBRARIES} ${MPI_CXX_LIBRARIES})
target_link_libraries(test_NCDFScaling ${MYLIBRARIES})
target_link_libraries(test_blocks  ${MYLIBRARIES})
target_link_libraries(test_poq cppunit ambitious)
target_link_libraries(test_timer cppunit ambitious)
target_link_libraries(test_bidirmap cppunit ambitious)
target_link_libraries(test_qtool cppunit ambitious)
target_link_libraries(test_po_tool cppunit ambitious ${Boost_LIBRARIES})
target_link_libraries(test_meshes ${MYLIBRARIES})
target_link_libraries(test_all ${MYLIBRARIES})
