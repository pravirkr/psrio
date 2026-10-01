include(FetchContent)

# Catch2: prefer system package (apt/brew) when available; pin 3.16.0 otherwise.
find_package(
  Catch2
  3
  CONFIG
  QUIET
)
if(Catch2_FOUND)
  message(STATUS "Using system Catch2 ${Catch2_VERSION}")
else()
  if(NOT DEFINED PSRIO_CATCH2_VERSION)
    set(PSRIO_CATCH2_VERSION "3.8.0")
  endif()
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v${PSRIO_CATCH2_VERSION}
    GIT_SHALLOW TRUE
  )
  set(CATCH_INSTALL_DOCS
      OFF
      CACHE BOOL "" FORCE
  )
  set(CATCH_INSTALL_EXTRAS
      OFF
      CACHE BOOL "" FORCE
  )
  FetchContent_MakeAvailable(Catch2)
endif()

if(NOT TARGET Catch2::Catch2WithMain)
  message(
    FATAL_ERROR "Catch2 is required for PSRIO_BUILD_TESTS=ON but was not found or downloaded."
  )
endif()

list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)

set(TARGET_TESTS ${PROJECT_NAME}_tests)

file(
  GLOB
  TEST_SOURCES
  CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

add_executable(${TARGET_TESTS} ${TEST_SOURCES})

target_link_libraries(
  ${TARGET_TESTS} PRIVATE ${PROJECT_NAME}::${PROJECT_NAME} Catch2::Catch2WithMain
                          psrio_developer_flags
)

target_compile_definitions(
  ${TARGET_TESTS} PRIVATE "PSRIO_TEST_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/data\""
)

include(Catch)
catch_discover_tests(${TARGET_TESTS})
