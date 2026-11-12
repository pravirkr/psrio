cmake_minimum_required(VERSION 3.21...3.30)

project(
  psrio
  VERSION 0.5.0
  LANGUAGES CXX
  DESCRIPTION "High-Performance Header-Only Pulsar & FRB I/O Library"
)

# Compiler version requirements
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11.0)
    message(
      FATAL_ERROR "GCC 11.0+ is required for C++20 support. Found: ${CMAKE_CXX_COMPILER_VERSION}."
    )
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0)
    message(
      FATAL_ERROR "Clang 15.0+ is required for C++20 support. Found: ${CMAKE_CXX_COMPILER_VERSION}."
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
option(PSRIO_BUILD_APPLICATIONS "Build psrio CLI applications" ${PSRIO_IS_TOP_LEVEL})
option(PSRIO_BUILD_DOCS "Build documentation" OFF)
option(PSRIO_ENABLE_COVERAGE "Enable code coverage for tests" OFF)
option(PSRIO_WITH_HDF5 "Enable the opt-in FBH5 reader" OFF)
option(PSRIO_WITH_PSRDADA "Enable the opt-in PSRDADA ring reader" OFF)
option(PSRIO_INSTALL "Generate the install rules" ${PSRIO_IS_TOP_LEVEL})

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

include(cmake/CPM.cmake)
if(DEFINED ENV{CPM_SOURCE_CACHE})
  set(CPM_SOURCE_CACHE
      $ENV{CPM_SOURCE_CACHE}
      CACHE PATH "Directory to cache downloaded CPM sources"
  )
endif()

# fmt (header-only) formats diagnostic and metadata strings. libstdc++ before
# GCC 13 has no std::format. A parent that already defines fmt::fmt-header-only
# (fmt 12.1 or newer) shares that target; define it before adding psrio.
# Otherwise fmt 12.1.0 is fetched. psrio does not install fmt, except with the
# applications: dmt installs spdlog unconditionally, and spdlog's export set
# needs fmt in one too.
if(PSRIO_BUILD_APPLICATIONS)
  set(PSRIO_FMT_INSTALL ON)
else()
  set(PSRIO_FMT_INSTALL OFF)
endif()
if(NOT TARGET fmt::fmt-header-only)
  CPMAddPackage(
    NAME fmt
    VERSION 12.1.0
    URL https://github.com/fmtlib/fmt/archive/refs/tags/12.1.0.tar.gz
    URL_HASH SHA256=ea7de4299689e12b6dddd392f9896f08fb0777ac7168897a244a6d6085043fea
    OPTIONS "FMT_INSTALL ${PSRIO_FMT_INSTALL}"
    EXCLUDE_FROM_ALL YES
    SYSTEM YES
  )
endif()
if(NOT TARGET fmt::fmt-header-only)
  message(FATAL_ERROR "fmt::fmt-header-only target not available after CPM setup.")
endif()
# Only the build interface names the in-tree target. Naming it from the install
# interface makes install(EXPORT) require fmt to be exported. The installed
# config appends fmt::fmt-header-only after find_dependency(fmt).
target_link_libraries(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:fmt::fmt-header-only>)

if(PSRIO_WITH_HDF5)
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
  add_library(psrio_psrdada INTERFACE)
  add_library(psrio::psrdada ALIAS psrio_psrdada)
  set_target_properties(psrio_psrdada PROPERTIES EXPORT_NAME psrdada)
  target_compile_features(psrio_psrdada INTERFACE cxx_std_20)
  if(PSRDADA_INCLUDE_DIR AND PSRDADA_LIBRARY)
    target_include_directories(psrio_psrdada INTERFACE ${PSRDADA_INCLUDE_DIR})
    target_link_libraries(psrio_psrdada INTERFACE ${PROJECT_NAME} ${PSRDADA_LIBRARY})
  else()
    if(NOT TARGET psrdada)
      enable_language(C)
      set(CMAKE_DISABLE_FIND_PACKAGE_CUDAToolkit TRUE)
      set(CMAKE_DISABLE_FIND_PACKAGE_CUDA TRUE)
      set(BUILD_TESTING OFF)
      CPMAddPackage(
        NAME psrdada
        GIT_REPOSITORY https://git.code.sf.net/p/psrdada/code
        GIT_TAG 151198ea3ec4a3a237f51e6217a5a2d7ff39194e
        EXCLUDE_FROM_ALL YES
        SYSTEM YES
      )
      unset(CMAKE_DISABLE_FIND_PACKAGE_CUDAToolkit)
      unset(CMAKE_DISABLE_FIND_PACKAGE_CUDA)
      unset(BUILD_TESTING)
    endif()
    if(NOT TARGET psrdada)
      message(
        FATAL_ERROR
          "PSRIO_WITH_PSRDADA=ON but PSRDADA was not found (need dada_hdu.h and libpsrdada)."
      )
    endif()
    target_link_libraries(psrio_psrdada INTERFACE ${PROJECT_NAME} $<BUILD_INTERFACE:psrdada>)
  endif()
endif()

# ------------------------------------------------------------------------------
# Testing, Applications & Developer Builds
# ------------------------------------------------------------------------------
if(PSRIO_BUILD_TESTS OR PSRIO_BUILD_APPLICATIONS)
  # Developer warning & optimization flags - isolated to test and CLI targets
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
    include(CheckCXXCompilerFlag)
    # Probe the positive flag. Clang accepts an unknown -Wno-* option silently,
    # so that form cannot tell Clang 15 from a compiler that has the group.
    check_cxx_compiler_flag("-Wc2y-extensions" PSRIO_HAS_C2Y_EXTENSIONS_WARNING)
    if(PSRIO_HAS_C2Y_EXTENSIONS_WARNING)
      target_compile_options(psrio_developer_flags INTERFACE -Wno-c2y-extensions)
    endif()
  endif()

  if(PSRIO_ENABLE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(psrio_developer_flags INTERFACE --coverage -O0 -g)
    target_link_options(psrio_developer_flags INTERFACE --coverage)
  endif()
endif()

# Applications come first so the tests can link the CLI command library.
if(PSRIO_BUILD_APPLICATIONS)
  add_subdirectory(applications)
endif()

if(PSRIO_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

# ------------------------------------------------------------------------------
# Installation & CMake Packaging (Header-Only / Arch-Independent)
# ------------------------------------------------------------------------------
if(PSRIO_INSTALL)
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
endif()

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
message(STATUS "--------------------------------------------------------")
message(STATUS "psrio ${PROJECT_VERSION} configured successfully")
message(STATUS "  - Type: Header-only library")
message(STATUS "  - Standard: C++20")
message(STATUS "  - Build Tests: ${PSRIO_BUILD_TESTS}")
message(STATUS "  - Build Applications: ${PSRIO_BUILD_APPLICATIONS}")
message(STATUS "  - FBH5 (HDF5): ${PSRIO_WITH_HDF5}")
message(STATUS "  - PSRDADA: ${PSRIO_WITH_PSRDADA}")
message(STATUS "--------------------------------------------------------")
