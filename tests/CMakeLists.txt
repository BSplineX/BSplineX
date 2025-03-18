include(BSplineXCatch2)
include(BSplineXJSON)

file(GLOB_RECURSE KNOTS_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/knots/test_*.cpp")
add_executable(test_knots ${KNOTS_TESTS})
target_link_libraries(test_knots PRIVATE BSplineX Catch2::Catch2WithMain)
catch_discover_tests(test_knots)

file(GLOB_RECURSE
  CONTROL_POINTS_TESTS
  "${CMAKE_CURRENT_SOURCE_DIR}/control_points/test_*.cpp"
)
add_executable(test_control_points ${CONTROL_POINTS_TESTS})
target_link_libraries(test_control_points PRIVATE
  BSplineX Catch2::Catch2WithMain
)
catch_discover_tests(test_control_points)

file(GLOB_RECURSE
  DEBOOR_TESTS
  "${CMAKE_CURRENT_SOURCE_DIR}/deboor/test_*.cpp"
)
add_executable(test_deboor ${DEBOOR_TESTS})
target_link_libraries(test_deboor PRIVATE
  BSplineX Catch2::Catch2WithMain
)
catch_discover_tests(test_deboor)

file(GLOB_RECURSE
  BSPLINE_TESTS
  "${CMAKE_CURRENT_SOURCE_DIR}/bspline/test_bspline.cpp"
)
add_executable(test_bspline ${BSPLINE_TESTS})
target_link_libraries(test_bspline PRIVATE
  BSplineX Catch2::Catch2WithMain
)
catch_discover_tests(test_bspline)

file(GLOB_RECURSE
  BSPLINE_TESTS_WITH_DATA
  "${CMAKE_CURRENT_SOURCE_DIR}/bspline/test_bspline_with_data.cpp"
)
add_executable(test_bspline_with_data ${BSPLINE_TESTS_WITH_DATA})
target_link_libraries(test_bspline_with_data PRIVATE
  BSplineX Catch2::Catch2WithMain nlohmann_json::nlohmann_json
)
catch_discover_tests(test_bspline_with_data)
