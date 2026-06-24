cmake_minimum_required(VERSION 3.16)
project(autoland_dev LANGUAGES CXX)

# -----------------------------------------------------------------------------
# Standalone autoland development environment.
# Every third-party dependency is pulled with FetchContent so the project builds
# with no system installs (network access is required on first configure).
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

# --- Eigen (header-only linear algebra) --------------------------------------
# Eigen's own CMake build is heavy; we only need the headers, so we populate the
# source and expose it as an INTERFACE target.
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
  GIT_SHALLOW TRUE
)
FetchContent_GetProperties(eigen)
if(NOT eigen_POPULATED)
  FetchContent_Populate(eigen)
endif()
add_library(eigen_headers INTERFACE)
target_include_directories(eigen_headers INTERFACE ${eigen_SOURCE_DIR})

# --- yaml-cpp (config parsing) -----------------------------------------------
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
set(YAML_CPP_INSTALL OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG 0.8.0
  GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(yaml-cpp)

# --- Catch2 (unit tests) -----------------------------------------------------
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.5.2
  GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(Catch2)

# -----------------------------------------------------------------------------
# Core library: every module lives here so tests, the app, and (later) a
# nonlinear plant or CBF-QP all link against the same single source of truth.
# -----------------------------------------------------------------------------
add_library(autoland_core
  src/aero_table.cpp
  src/config.cpp
  src/mixing.cpp
  src/dynamics.cpp
  src/trim.cpp
  src/linear_model.cpp
  src/controller.cpp
  src/cbf.cpp
  src/sim.cpp
)
target_include_directories(autoland_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(autoland_core PUBLIC eigen_headers yaml-cpp::yaml-cpp)
target_compile_options(autoland_core PRIVATE -Wall -Wextra)

# Path to the bundled data directory, so the app/tests can find configs by default.
target_compile_definitions(autoland_core PUBLIC
  AUTOLAND_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")

# -----------------------------------------------------------------------------
# Simulation application
# -----------------------------------------------------------------------------
add_executable(autoland_sim apps/autoland_sim.cpp)
target_link_libraries(autoland_sim PRIVATE autoland_core)

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------
enable_testing()
include(CTest)
include(Catch)

add_executable(unit_tests
  test/test_aero_table.cpp
  test/test_trim.cpp
  test/test_linear_model.cpp
)
target_link_libraries(unit_tests PRIVATE autoland_core Catch2::Catch2WithMain)
catch_discover_tests(unit_tests)
