cmake_minimum_required(VERSION 3.21...3.30)

project(
  psrio
  VERSION 0.3.0
  LANGUAGES CXX
  DESCRIPTION "High-Performance Header-Only Pulsar & FRB I/O Library"
)

# Compiler version requirements
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13.0)
    message(
      FATAL_ERROR "GCC 13.0+ is required for C++20 support. Found: ${CMAKE_CXX_COMPILER_VERSION}."
    )
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0)
    message(
      FATAL_ERROR "Clang 18.0+ is required for C++20 support. Found: ${CMAKE_CXX_COMPILER_VERSION}."
    )
  endif()
else()
  message(
    WARNING
      "Compiler '${CMAKE_CXX_COMPILER_ID}' is not officially tested. Proceed at your own risk."
  )
endif()

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(PSRIO_IS_TOP_LEVEL TRUE)
else()
  set(PSRIO_IS_TOP_LEVEL FALSE)
endif()

# Configuration options
option(PSRIO_BUILD_TESTS "Build test suite" ${PSRIO_IS_TOP_LEVEL})
option(PSRIO_BUILD_DOCS "Build documentation" OFF)
option(PSRIO_ENABLE_COVERAGE "Enable code coverage for tests" OFF)
option(PSRIO_WITH_HDF5 "Enable the opt-in FBH5 reader" OFF)
option(PSRIO_WITH_PSRDADA "Enable the opt-in PSRDADA ring reader" OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ------------------------------------------------------------------------------
# Core Interface Target Definition
# ------------------------------------------------------------------------------
include(GNUInstallDirs)

add_library(${PROJECT_NAME} INTERFACE)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_20)

target_include_directories(
  ${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

if(PSRIO_WITH_HDF5)
  include(cmake/CPM.cmake)
  if(DEFINED ENV{CPM_SOURCE_CACHE})
    set(CPM_SOURCE_CACHE
        $ENV{CPM_SOURCE_CACHE}
        CACHE PATH "Directory to cache downloaded CPM sources"
    )
  endif()

  CPMAddPackage(
    NAME HighFive
    VERSION 3.3.0
    URL https://github.com/highfive-devs/highfive/archive/refs/tags/v3.3.0.tar.gz
    OPTIONS "HIGHFIVE_FIND_HDF5 ON"
    EXCLUDE_FROM_ALL YES
    SYSTEM YES
  )

  if(NOT TARGET HighFive::HighFive)
    find_package(HighFive QUIET)
  endif()

  if(NOT TARGET HighFive::HighFive)
    message(FATAL_ERROR "HighFive::HighFive target not available after CPM setup.")
  endif()

  add_library(psrio_fbh5 INTERFACE)
  add_library(psrio::fbh5 ALIAS psrio_fbh5)
  set_target_properties(psrio_fbh5 PROPERTIES EXPORT_NAME fbh5)
  target_compile_features(psrio_fbh5 INTERFACE cxx_std_20)
  target_link_libraries(psrio_fbh5 INTERFACE ${PROJECT_NAME} HighFive::HighFive)
endif()

if(PSRIO_WITH_PSRDADA)
  find_path(PSRDADA_INCLUDE_DIR NAMES dada_hdu.h)
  find_library(PSRDADA_LIBRARY NAMES psrdada)
  if(NOT PSRDADA_INCLUDE_DIR OR NOT PSRDADA_LIBRARY)
    message(
      FATAL_ERROR
        "PSRIO_WITH_PSRDADA=ON but PSRDADA was not found (need dada_hdu.h and libpsrdada)."
    )
  endif()
  add_library(psrio_psrdada INTERFACE)
  add_library(psrio::psrdada ALIAS psrio_psrdada)
  set_target_properties(psrio_psrdada PROPERTIES EXPORT_NAME psrdada)
  target_compile_features(psrio_psrdada INTERFACE cxx_std_20)
  target_include_directories(psrio_psrdada INTERFACE ${PSRDADA_INCLUDE_DIR})
  target_link_libraries(psrio_psrdada INTERFACE ${PROJECT_NAME} ${PSRDADA_LIBRARY})
endif()

# ------------------------------------------------------------------------------
# Testing & Developer Builds
# ------------------------------------------------------------------------------
if(PSRIO_BUILD_TESTS)
  enable_testing()

  # Developer warning & optimization flags - isolated to test targets
  add_library(psrio_developer_flags INTERFACE)
  target_compile_options(
    psrio_developer_flags
    INTERFACE -Wall
              -Wextra
              -Wpedantic
              $<$<CONFIG:Debug>:-O0
              -g>
              $<$<CONFIG:Release>:-O3>
  )

  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(psrio_developer_flags INTERFACE -Wno-c2y-extensions)
  endif()

  if(PSRIO_ENABLE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(psrio_developer_flags INTERFACE --coverage -O0 -g)
    target_link_options(psrio_developer_flags INTERFACE --coverage)
  endif()

  add_subdirectory(tests)
endif()

# ------------------------------------------------------------------------------
# Installation & CMake Packaging (Header-Only / Arch-Independent)
# ------------------------------------------------------------------------------
include(CMakePackageConfigHelpers)

set(PSRIO_CMAKE_PACKAGE_DIR "${CMAKE_INSTALL_DATADIR}/cmake/${PROJECT_NAME}")
set(PSRIO_INSTALL_TARGETS ${PROJECT_NAME})
if(TARGET psrio_fbh5)
  list(APPEND PSRIO_INSTALL_TARGETS psrio_fbh5)
endif()
if(TARGET psrio_psrdada)
  list(APPEND PSRIO_INSTALL_TARGETS psrio_psrdada)
endif()
install(
  TARGETS ${PSRIO_INSTALL_TARGETS}
  EXPORT psrioTargets
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(
  EXPORT psrioTargets
  FILE psrioTargets.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${PSRIO_CMAKE_PACKAGE_DIR}
)

configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/psrioConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/psrioConfig.cmake
  INSTALL_DESTINATION ${PSRIO_CMAKE_PACKAGE_DIR}
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/psrioConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion ARCH_INDEPENDENT
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/psrioConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/psrioConfigVersion.cmake
        DESTINATION ${PSRIO_CMAKE_PACKAGE_DIR}
)

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
message(STATUS "--------------------------------------------------------")
message(STATUS "psrio ${PROJECT_VERSION} configured successfully")
message(STATUS "  - Type: Header-only library")
message(STATUS "  - Standard: C++20")
message(STATUS "  - Build Tests: ${PSRIO_BUILD_TESTS}")
message(STATUS "  - FBH5 (HDF5): ${PSRIO_WITH_HDF5}")
message(STATUS "  - PSRDADA: ${PSRIO_WITH_PSRDADA}")
message(STATUS "--------------------------------------------------------")
