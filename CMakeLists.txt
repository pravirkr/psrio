cmake_minimum_required(VERSION 3.21...3.30)

project(
  psrio
  VERSION 0.1.0
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
option(PSRIO_BUILD_TESTS "Build test suite" ${PROJECT_IS_TOP_LEVEL})
option(PSRIO_BUILD_DOCS "Build documentation" OFF)
option(PSRIO_ENABLE_COVERAGE "Enable code coverage for tests" OFF)

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
install(
  TARGETS ${PROJECT_NAME}
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
message(STATUS "--------------------------------------------------------")
