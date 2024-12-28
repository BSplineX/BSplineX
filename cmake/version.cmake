# This file will fetch the version from git tags and use it in CMake. Adapted
# from
# https://github.com/behnamasadi/cpp_tutorials/blob/980e82c92d89bb9096eff63cd5e23a50b4bd3846/docs/getting_version_from_git_in_CMake.md

find_package(Git QUIET)

set(GIT_VERSION_AVAILABLE TRUE)
if(NOT GIT_FOUND)
  set(GIT_VERSION_AVAILABLE FALSE)
else()
  # Get the latest tag
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE BSPLINEX_LATEST_RELEASE
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if("${BSPLINEX_LATEST_RELEASE}" STREQUAL "")
    set(GIT_VERSION_AVAILABLE FALSE)
  else()
    # Get the number of commits since the last tag
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-list ${BSPLINEX_LATEST_RELEASE}..HEAD --count
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE GIT_COMMITS_SINCE_TAG
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    if("${GIT_COMMITS_SINCE_TAG}" STREQUAL "0")
      set(BSPLINEX_VERSION "${BSPLINEX_LATEST_RELEASE}")
      message(STATUS "BSplineX: Release version ${BSPLINEX_VERSION}")
    else()
      set(BSPLINEX_VERSION "${BSPLINEX_LATEST_RELEASE}.${GIT_COMMITS_SINCE_TAG}")
      message(STATUS "BSplineX: Development version ${BSPLINEX_VERSION}")
    endif()
  endif()
endif()

if(NOT GIT_VERSION_AVAILABLE)
  set(BSPLINEX_VERSION "0.0.0")
  message(STATUS "BSplineX: Could not determine version from git tags. Using ${BSPLINEX_VERSION}")
endif()

set(BSPLINEX_VERSION_FULL
  "${BSPLINEX_VERSION}"
  CACHE INTERNAL "BSplineX version")
