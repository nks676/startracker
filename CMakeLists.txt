cmake_minimum_required(VERSION 3.14)
project(startracker VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# We will need to include src directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# Define the executable
add_executable(startracker
    src/main.cpp
    src/image_processing.cpp
    src/catalog.cpp
    src/identification.cpp
    src/estimation.cpp
)

# STB Image for reading test images
# Eigen for matrix math (we will download it or require it to be installed)
# For now, let's just make it compile with dummy files.
