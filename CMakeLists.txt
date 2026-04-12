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
    src/identification.cpp
    src/estimation.cpp
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
)
target_link_libraries(startracker_tests PRIVATE startracker_core GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(startracker_tests)
