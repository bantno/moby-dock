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

# --- autodiff (exact higher-order forward-mode derivatives for the HOCBFs) ----
# Header-only C++17, used by the augmented-longitudinal CBF to compute Lie
# derivatives exactly (no finite differences). autodiff's own CMake does a
# find_package(Eigen3), which we can't satisfy with our header-only Eigen, so we
# only POPULATE the source and expose its include dir as an INTERFACE target
# (same approach as Eigen above). Headers live under ${autodiff_SOURCE_DIR}.
FetchContent_Declare(
  autodiff
  GIT_REPOSITORY https://github.com/autodiff/autodiff.git
  GIT_TAG v1.1.2
  GIT_SHALLOW TRUE
)
FetchContent_GetProperties(autodiff)
if(NOT autodiff_POPULATED)
  FetchContent_Populate(autodiff)
endif()
add_library(autodiff_headers INTERFACE)
target_include_directories(autodiff_headers INTERFACE ${autodiff_SOURCE_DIR})
target_link_libraries(autodiff_headers INTERFACE eigen_headers)

# --- OSQP (quadratic-program solver backing the CBF safety filter) -----------
# Linked privately into autoland_core; only qp_solver.cpp includes <osqp.h>, so
# OSQP is an implementation detail and does not leak into the public headers.
set(OSQP_BUILD_DEMO_EXE OFF CACHE BOOL "" FORCE)
set(UNITTESTS OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
  osqp
  GIT_REPOSITORY https://github.com/osqp/osqp.git
  GIT_TAG v0.6.3
  GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(osqp)

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
  src/nominal_lon_glideslope.cpp
  src/qp_solver.cpp
  src/cbf.cpp
  src/sim.cpp
  src/lon_augmented.cpp
  src/lon_cbf_filter.cpp
  src/lon_sim.cpp
)
target_include_directories(autoland_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# autodiff is used in public headers (lie_taylor.hpp), so it must be PUBLIC.
target_link_libraries(autoland_core PUBLIC eigen_headers yaml-cpp::yaml-cpp autodiff_headers)
target_link_libraries(autoland_core PRIVATE osqpstatic)
# Warning flags are GCC/Clang spelling; MSVC uses /W and doesn't define M_PI
# from <cmath> unless _USE_MATH_DEFINES is set. Keep the project buildable on
# both toolchains without weakening the GCC/Clang build.
if(MSVC)
  target_compile_options(autoland_core PRIVATE /W3)
  target_compile_definitions(autoland_core PUBLIC _USE_MATH_DEFINES)
else()
  target_compile_options(autoland_core PRIVATE -Wall -Wextra)
endif()

# Path to the bundled data directory, so the app/tests can find configs by default.
target_compile_definitions(autoland_core PUBLIC
  AUTOLAND_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")

# -----------------------------------------------------------------------------
# Simulation application
# -----------------------------------------------------------------------------
add_executable(autoland_sim apps/autoland_sim.cpp)
target_link_libraries(autoland_sim PRIVATE autoland_core)

add_executable(lon_autoland_sim apps/lon_autoland_sim.cpp)
target_link_libraries(lon_autoland_sim PRIVATE autoland_core)

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------
enable_testing()
include(CTest)
include(Catch)

add_executable(unit_tests
  test/test_aero_table.cpp
  test/test_dynamics.cpp
  test/test_trim.cpp
  test/test_linear_model.cpp
  test/test_cbf.cpp
  test/test_lon_cbf.cpp
  test/test_lon_sim.cpp
  test/test_lie_derivatives.cpp
  test/test_stall_model.cpp
  test/test_wind_gust.cpp
  test/test_stall_recovery.cpp
)
target_link_libraries(unit_tests PRIVATE autoland_core Catch2::Catch2WithMain)
catch_discover_tests(unit_tests)
