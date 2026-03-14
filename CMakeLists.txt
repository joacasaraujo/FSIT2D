cmake_minimum_required(VERSION 3.10)

project(FSIT)

# Set the path for the Eigen library
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

# Include the Eigen directory
include_directories(${EIGEN3_DIR})

set(ENTRY_FILE "scenarios/validation.cpp" CACHE STRING "scenarios/Entry point .cpp file")

set(FSIT_SOURCES
    ${ENTRY_FILE}
    lbm/Cell.cpp
    dem/Body.cpp
    dem/Interaction.cpp
    engine/Engine.cpp
    property/Shape.cpp
    scene/Scene.cpp
    scene/Output.cpp)

# Define your executable and its source files
add_executable(FSIT ${FSIT_SOURCES})

# Set the include directories for the FSIT target
target_include_directories(FSIT PUBLIC ${PROJECT_SOURCE_DIR})

target_link_libraries (FSIT Eigen3::Eigen)
