cmake_minimum_required(VERSION 3.2)
project(D2K3Pi VERSION 0.1.0)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ROOT enforces at least C++14
# Something else, possibly Boost, enforces C++17
set(CMAKE_CXX_STANDARD 17)

# ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Minuit2 RooFit RooStats RooFitCore)
include(${ROOT_USE_FILE})

# Boost
# Need at least version 1.66.0 for numerical integration
find_package( Boost 1.66.0 REQUIRED COMPONENTS system filesystem)
include_directories(${Boost_INCLUDE_DIRS})

# Some standard flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wformat-security -Werror -O0 -g")

# Debug
set(CMAKE_BUILD_TYPE Debug)

# Common functionality like plotting helpers + stuff
add_subdirectory("${CMAKE_SOURCE_DIR}/common")

# Add our fitter library
add_subdirectory("${CMAKE_SOURCE_DIR}/charmFitter")

# Add an efficiency library
add_subdirectory("${CMAKE_SOURCE_DIR}/efficiency")

# Library for doing sWeighting
add_subdirectory("${CMAKE_SOURCE_DIR}/lhcbMonteCarlo/sWeighting")

# UT
enable_testing()
add_subdirectory("${CMAKE_SOURCE_DIR}/test/ut/")

# IT
add_subdirectory("${CMAKE_SOURCE_DIR}/test/it/")
