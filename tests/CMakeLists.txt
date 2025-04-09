include(BSplineXCatch2)
include(BSplineXJSON)

file(GLOB_RECURSE KNOTS_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/knots/test_*.cpp")
add_executable(test_knots ${KNOTS_TESTS})
target_include_directories(test_knots PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_knots PRIVATE BSplineX Catch2::Catch2WithMain)
catch_discover_tests(test_knots)

file(GLOB_RECURSE
  CONTROL_POINTS_TESTS
  "${CMAKE_CURRENT_SOURCE_DIR}/control_points/test_*.cpp"
)
add_executable(test_control_points ${CONTROL_POINTS_TESTS})
target_include_directories(test_control_points PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_control_points PRIVATE
  BSplineX Catch2::Catch2WithMain
)
catch_discover_tests(test_control_points)

file(GLOB_RECURSE
  BSPLINE_TESTS_WITH_DATA
  "${CMAKE_CURRENT_SOURCE_DIR}/bspline/test_bspline_with_data.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/bspline/matchers.hpp"
)
add_executable(test_bspline_with_data ${BSPLINE_TESTS_WITH_DATA})
target_include_directories(test_bspline_with_data PRIVATE "${CMAKE_SOURCE_DIR}/tests/")
target_link_libraries(test_bspline_with_data PRIVATE
  BSplineX Catch2::Catch2WithMain nlohmann_json::nlohmann_json
)
target_compile_definitions(test_bspline_with_data PRIVATE
  TEST_RESOURCES_DIR="${CMAKE_SOURCE_DIR}/tests/bspline/data"
)
catch_discover_tests(test_bspline_with_data)
