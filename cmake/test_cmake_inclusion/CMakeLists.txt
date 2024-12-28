cmake_minimum_required(VERSION 3.14)

project(test)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

option(ADD_SUBDIRECTORY_METHOD OFF)
option(FETCH_CONTENT_METHOD ON)
option(FIND_PACKAGE_METHOD OFF)

set(CMAKE_BUILD_TYPE Release)

add_executable(test test.cpp)

if(ADD_SUBDIRECTORY_METHOD)
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/third_party/BSplineX")
endif()

if(FIND_PACKAGE_METHOD)
  list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
  find_package(
    BSplineX
    0.0.0
    NO_MODULE
  )
endif()

if(FETCH_CONTENT_METHOD)
  include(FetchContent)

  set(FETCHCONTENT_BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
  fetchcontent_declare(
      BSplineX
      GIT_REPOSITORY https://github.com/BSplineX/BSplineX.git
      GIT_TAG        cmake_tests
  )

  fetchcontent_makeavailable(BSplineX)
endif()

target_link_libraries(test PRIVATE BSplineX::BSplineX)

