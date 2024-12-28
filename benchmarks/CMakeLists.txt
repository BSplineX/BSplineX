include(catch)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

# Benchmarks

file(GLOB_RECURSE
  BSPLINE_BENCHS
  "${CMAKE_CURRENT_SOURCE_DIR}/bspline/bench_*.cpp"
)
add_executable(bench_bspline ${BSPLINE_BENCHS})

target_link_libraries(bench_bspline PRIVATE
  BSplineX Catch2::Catch2WithMain
)
catch_discover_tests(bench_bspline)

