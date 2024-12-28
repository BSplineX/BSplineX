function(parse_version_string VERSION_STRING)
  # Remove 'v' prefix if present
  string(REGEX REPLACE "^v" "" VERSION_STRING "${VERSION_STRING}")

  # pre-release tag example: 0.0.1-alpha.1
  string(REGEX MATCH "^([0-9]+\.[0-9]+\.[0-9]+)(-[a-zA-Z0-9]+)?$" _ "${VERSION_STRING}")
  set(BASE_VERSION "${CMAKE_MATCH_1}")
  set(PRE_RELEASE "${CMAKE_MATCH_2}")

  if(PRE_RELEASE)
    string(SUBSTRING "${PRE_RELEASE}" 1 -1 PRE_RELEASE) # Remove the leading "-"
  endif()

  set(BSPLINEX_VERSION ${BASE_VERSION} CACHE INTERNAL "")
  set(BSPLINEX_VERSION_PRERELEASE ${PRE_RELEASE} CACHE INTERNAL "")
endfunction()

function(get_version_from_git)
  # Try to get version from git first
  find_package(Git QUIET)
  if(GIT_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse --is-inside-work-tree
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE IS_GIT_REPO
      OUTPUT_QUIET
      ERROR_QUIET
    )

    if(IS_GIT_REPO EQUAL 0)
      execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_TAG
        ERROR_VARIABLE GIT_ERROR
        RESULT_VARIABLE GIT_RESULT
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )

      if(GIT_RESULT EQUAL 0)
        parse_version_string("${GIT_TAG}")
        return()
      endif()
    endif()
  endif()

  # Fallback: Try to read version from version file
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION")
    file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" FILE_VERSION)
    string(STRIP "${FILE_VERSION}" FILE_VERSION)
    parse_version_string("${FILE_VERSION}")
    return()
  endif()

  # Final fallback: Use a default version
  set(BSPLINEX_VERSION "0.0.0" CACHE INTERNAL "")
  set(BSPLINEX_VERSION_PRERELEASE "" CACHE INTERNAL "")
endfunction()