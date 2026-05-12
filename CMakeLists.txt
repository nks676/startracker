cmake_minimum_required(VERSION 3.14)
project(startracker VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include src directory for headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# --- Core library (shared by main executable and tests) ---
add_library(startracker_core STATIC
    src/image_processing.cpp
    src/catalog.cpp
    src/camera_model.cpp
    src/wahba.cpp
    src/pattern_hash.cpp
    src/identification.cpp
    src/estimation.cpp
    src/tiff_reader.cpp
)

# --- Main executable ---
add_executable(startracker src/main.cpp)
target_link_libraries(startracker PRIVATE startracker_core)

# --- Testing (GoogleTest via FetchContent) ---
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.14.0
)
FetchContent_MakeAvailable(googletest)
enable_testing()

add_executable(startracker_tests
    tests/test_image_processing.cpp
    tests/test_catalog.cpp
    tests/test_identification.cpp
    tests/test_estimation.cpp
    tests/test_camera_model.cpp
)
target_link_libraries(startracker_tests PRIVATE startracker_core GTest::gtest_main)

include(GoogleTest)
# DISCOVERY_MODE PRE_TEST: enumerate test cases when `ctest` runs rather than
# at build time. POST_BUILD (the default) invokes the freshly-linked test
# binary as a build step, which fails under cross-compile because the host
# can't execute the aarch64 ELF — surfaces as
# "GoogleTestAddTests: Unknown system error -8" with ENOEXEC. The host CI
# lane sees no behavioural change (ctest still discovers + runs everything);
# the build-arm64-cross lane now finishes the compile-and-link signal it
# was meant to deliver. CMake documents PRE_TEST as the cross-compile fix.
gtest_discover_tests(startracker_tests DISCOVERY_MODE PRE_TEST)
