set(CATCH2_REQUIRED_VERSION 3.7.1)
cmake_policy(SET CMP0135 NEW)

list(APPEND CMAKE_PATH_PREFIX "${BSPLINEX_THIRD_PARTY_DIR}")
find_package(
  Catch2
  ${CATCH2_REQUIRED_VERSION}
  NO_MODULE
  QUIET
)

if(NOT TARGET Catch2::Catch2WithMain)
  # We have to do this manually since the install process of
  # Catch2 does not make the target directly discoverable by find_package
  if(EXISTS "${BSPLINEX_THIRD_PARTY_DIR}/catch2-src")
    message(STATUS
      "BSplineX: "
      "Found Catch2 installed in ${BSPLINEX_THIRD_PARTY_DIR}/catch2-src"
    )
    add_subdirectory(
      "${BSPLINEX_THIRD_PARTY_DIR}/catch2-src"
      "${BSPLINEX_THIRD_PARTY_DIR}/catch2-build"
    )
  else()
    message(STATUS "BSplineX: "
      "Did not find Catch2 ${CATCH2_REQUIRED_VERSION} installed, "
      "downloading to ${BSPLINEX_THIRD_PARTY_DIR}"
    )
    include(FetchContent)

    set(FETCHCONTENT_BASE_DIR "${BSPLINEX_THIRD_PARTY_DIR}")
    fetchcontent_declare(
      Catch2
      URL "https://github.com/catchorg/Catch2/archive/refs/tags/v${CATCH2_REQUIRED_VERSION}.tar.gz"
    )

    fetchcontent_makeavailable(Catch2)
  endif()
endif()

get_target_property(CATCH2_SOURCE_DIR
  Catch2::Catch2
  SOURCE_DIR
)

if(NOT DEFINED CATCH2_SOURCE_DIR)
  message(FATAL_ERROR "BSplineX: Could not find Catch2 source directory")
endif()

list(APPEND CMAKE_MODULE_PATH ${CATCH2_SOURCE_DIR}/../extras)
include(CTest)
include(Catch)