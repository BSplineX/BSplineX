cmake_minimum_required(VERSION 3.14)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(BSplineXVersion)
get_version_from_git()

project(BSplineX
  VERSION ${BSPLINEX_VERSION}
  DESCRIPTION "A header-only library to harness the power of BSplines"
  HOMEPAGE_URL "https://github.com/BSplineX/BSplineX"
  LANGUAGES CXX
)

if(NOT BSPLINEX_VERSION_PRERELEASE)
  set(BSPLINEX_FULL_VERSION "${BSPLINEX_VERSION}")
  write_version_file()
else()
  set(BSPLINEX_FULL_VERSION "${BSPLINEX_VERSION}-${BSPLINEX_VERSION_PRERELEASE}")
endif()
message(STATUS "BSplineX: version ${BSPLINEX_FULL_VERSION}")

set(BSPLINEX_THIRD_PARTY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

if(NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib")
  message(STATUS
    "BSplineX: "
    "CMAKE_INSTALL_PREFIX not specified, "
    "defaulting to ${CMAKE_CURRENT_SOURCE_DIR}"
  )
endif()

# Dependencies
include(BSplineXEigen3)

# Library definition
add_library(BSplineX INTERFACE)
add_library(BSplineX::BSplineX ALIAS BSplineX)

target_link_libraries(BSplineX INTERFACE Eigen3::Eigen)

target_include_directories(BSplineX INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# INFO:
# Ensure this options are only available to developers.
# In theory we could leave them available, but the call to
# `target_sources` breaks `fetchcontent` compatibility due
# to absolute paths being added to `INTERFACE_SOURCES`.
# I tried solving it, but it seems to be poorly documented,
# supported, etc.
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
    message(STATUS
      "BSplineX: "
      "CMAKE_BUILD_TYPE not specified, defaulting to ${CMAKE_BUILD_TYPE}"
    )
  endif()

  option(BUILD_TESTS "Build tests" OFF)
  option(BUILD_BENCHMARKS "Build benchmarks" OFF)
  option(BUILD_EXAMPLES "Build examples" OFF)
  option(BSPLINEX_DEBUG_LOG_CALL "Debug class contructors/destructors" OFF)

  file(GLOB_RECURSE HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
  foreach(HEADER_FILE IN LISTS HEADER_FILES)
    file(
      RELATIVE_PATH
      RELATIVE_PATH_HEADER_FILE
      ${CMAKE_CURRENT_SOURCE_DIR} ${HEADER_FILE}
    )
    list(APPEND RELATIVE_HEADER_FILES ${RELATIVE_PATH_HEADER_FILE})
  endforeach()
  target_sources(BSplineX INTERFACE $<BUILD_INTERFACE:${RELATIVE_HEADER_FILES}>)

  if(BSPLINEX_DEBUG_LOG_CALL)
    target_compile_definitions(BSplineX INTERFACE -DBSPLINEX_DEBUG_LOG_CALL)
  endif()

  if(BUILD_TESTS)
    string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)
  endif()

  if(BUILD_BENCHMARKS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/benchmarks)
  endif()

  if(BUILD_EXAMPLES)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)
  endif()
endif()

if(MSVC)
  target_compile_options(BSplineX INTERFACE /W4 /WX)
else()
  target_compile_options(BSplineX INTERFACE -Wall -Wextra -pedantic -Werror)
endif()

# Installation
install(TARGETS BSplineX EXPORT BSplineXTargets)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION include)

include(CMakePackageConfigHelpers)

export(EXPORT BSplineXTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/BSplineXTargets.cmake
  NAMESPACE BSplineX::
)

install(EXPORT BSplineXTargets
  FILE BSplineXTargets.cmake
  NAMESPACE BSplineX::
  DESTINATION lib/cmake/BSplineX
)

write_basic_package_version_file(
  BSplineXConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/BSplineXConfig.cmake
  INSTALL_DESTINATION lib/cmake/BSplineX
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/BSplineXConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/BSplineXConfigVersion.cmake
  DESTINATION lib/cmake/BSplineX
)
