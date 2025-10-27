cmake_minimum_required(VERSION 3.15)
project(ChemPlasKin VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(Cantera REQUIRED)
find_package(PkgConfig REQUIRED)

# Find OpenMP if available
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_compile_options(ChemPlasKin PRIVATE OpenMP::OpenMP_CXX)
endif()

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Collect source files
file(GLOB_RECURSE SOURCES "src/*.cpp")
file(GLOB_RECURSE HEADERS "include/*.h")

# Create executable
add_executable(ChemPlasKin ${SOURCES} ${HEADERS})

# Link libraries
target_link_libraries(ChemPlasKin PRIVATE cantera_shared)

# Add compiler-specific options
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(ChemPlasKin PRIVATE -Wall -Wextra -O3)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(ChemPlasKin PRIVATE -Wall -Wextra -O3)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(ChemPlasKin PRIVATE /W3 /O2)
endif()

# Installation
install(TARGETS ChemPlasKin DESTINATION bin)
install(FILES ${HEADERS} DESTINATION include/ChemPlasKin)

# Testing
enable_testing()
add_subdirectory(tests OPTIONAL)
