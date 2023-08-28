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
        LANGUAGES  C CXX CUDA
        )

set(CMAKE_CUDA_ARCHITECTURES 80)
include_directories("${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")

include(CTest)
enable_testing()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

#
# Prevent building in the source directory
#

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n")
endif()
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

#find_package(ColPack REQUIRED)
find_package(OpenMP REQUIRED)

# If CUDA enabled, find nvcc
# if(ENABLE_CUDA)
#   find_package(CUDA REQUIRED)
#   set(CMAKE_CXX_COMPILER ${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc)
#   set(CMAKE_CUDA_COMPILER ${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc)
#   set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
#   set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 -Xcompiler -fopenmp")
#   # Set CUDA_ARCHITECTURES to 80
#   set(CMAKE_CUDA_ARCHITECTURES 80)
#   include_directories(${CUDA_INCLUDE_DIRS})
#   link_directories(${CUDA_LIBRARY_DIRS})
# endif()


# The executable code is here
add_subdirectory(apps)

# The compiled library code is here
add_subdirectory(src)

