cmake_minimum_required(VERSION 3.1.0)

if(NOT CMAKE_CXX_COMPILER)
  if(APPLE)
    set(CMAKE_CXX_COMPILER "g++-12")
  elseif(UNIX)
    set(CMAKE_CXX_COMPILER "g++")
  endif()
endif()

project( CliqePart VERSION 0.1.0
        DESCRIPTION "Cliqe Partition for quantum chemistry"
        LANGUAGES  C CXX
        )

# If CUDA detected, set ENABLE_CUDA flag
find_package(CUDA)
if(CUDA_FOUND)
  set(ENABLE_CUDA ON)
  # Print CUDA found
  message(STATUS "CUDA found, enabling CUDA compilation")
else()
  set(ENABLE_CUDA OFF)
  message(STATUS "CUDA not found, disabling CUDA compilation")
endif()

# If CUDA enabled, add as language
if(ENABLE_CUDA)
  enable_language(CUDA)
  set(CMAKE_CUDA_ARCHITECTURES 80)
  include_directories("${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
endif()

include(CTest)
enable_testing()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#
# Prevent building in the source directory
#

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n")
endif()
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

#find_package(ColPack REQUIRED)
find_package(OpenMP REQUIRED)


# The executable code is here
add_subdirectory(apps)

# The compiled library code is here
add_subdirectory(src)

