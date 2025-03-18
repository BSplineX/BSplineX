set(JSON_REQUIRED_VERSION 3.11.3)
cmake_policy(SET CMP0135 NEW)

list(APPEND CMAKE_PREFIX_PATH "${BSPLINEX_THIRD_PARTY_DIR}")
find_package(
  nlohmann_json
  ${JSON_REQUIRED_VERSION}
  NO_MODULE
  QUIET
)

if(NOT TARGET nlohmann_json::nlohmann_json)
  if(EXISTS "${BSPLINEX_THIRD_PARTY_DIR}/json-src")
    message(STATUS
      "BSplineX: "
      "Found JSON installed in ${BSPLINEX_THIRD_PARTY_DIR}/json-src"
    )

    # cmake-lint: disable=C0103
    set(JSON_BuildTests OFF CACHE INTERNAL "")
    add_subdirectory(
      "${BSPLINEX_THIRD_PARTY_DIR}/json-src"
      "${BSPLINEX_THIRD_PARTY_DIR}/json-build"
    )
  else()

    message(STATUS
      "BSplineX: "
      "Did not find JSON ${JSON_REQUIRED_VERSION} installed, "
      "downloading to ${BSPLINEX_THIRD_PARTY_DIR}"
    )
    include(FetchContent)

    set(FETCHCONTENT_BASE_DIR "${BSPLINEX_THIRD_PARTY_DIR}")
    fetchcontent_declare(
      json
      URL  "https://github.com/nlohmann/json/releases/download/v${JSON_REQUIRED_VERSION}/json.tar.xz"
    )

    fetchcontent_makeavailable(json)
  endif()

else()
  get_target_property(JSON_INCLUDE_DIRS
    nlohmann_json::nlohmann_json
    INTERFACE_INCLUDE_DIRECTORIES
  )
  message(STATUS "BSplineX: Found JSON installed in ${JSON_INCLUDE_DIRS}")
endif()


