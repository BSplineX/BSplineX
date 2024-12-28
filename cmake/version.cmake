# This file will fetch the version from git tags and use it in CMake. Adapted
# from
# https://github.com/behnamasadi/cpp_tutorials/blob/980e82c92d89bb9096eff63cd5e23a50b4bd3846/docs/getting_version_from_git_in_CMake.md

find_package(Git QUIET)

if(NOT GIT_FOUND)
  set(BSPLINEX_VERSION "v0.0.0")
  set(BSPLINEX_LATEST_RELEASE ${BSPLINEX_VERSION})
  message(
    STATUS "Git not found, setting BSplineX version to ${BSPLINEX_VERSION}")
else()
  # Get the latest tag
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE BSPLINEX_LATEST_RELEASE
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Get the number of commits since the last tag
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-list ${BSPLINEX_LATEST_RELEASE}..HEAD --count
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMITS_SINCE_TAG
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if("${GIT_COMMITS_SINCE_TAG}" STREQUAL "0")
    set(BSPLINEX_VERSION "${BSPLINEX_LATEST_RELEASE}")
    message(STATUS "BSplineX release version: ${BSPLINEX_VERSION}")
  else()
    set(BSPLINEX_VERSION "${BSPLINEX_LATEST_RELEASE}.${GIT_COMMITS_SINCE_TAG}")
    message(STATUS "BSplineX development version: ${BSPLINEX_VERSION}")
  endif()
endif()

set(BSPLINEX_VERSION_FULL
  "${BSPLINEX_VERSION}"
  CACHE INTERNAL "BSplineX version")
