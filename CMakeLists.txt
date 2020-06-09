cmake_minimum_required(VERSION 3.2)
project(D2K3Pi VERSION 0.1.0)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ROOT enforces C++14
set(CMAKE_CXX_STANDARD 17)

# ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Minuit2)
include(${ROOT_USE_FILE})

# Boost
# Need at least version 1.66.0 for numerical integration
find_package( Boost 1.66.0 REQUIRED COMPONENTS system filesystem)
include_directories(${Boost_INCLUDE_DIRS})

# Some standard flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wformat-security -Werror -Og -g")

# Add our core library
# This defines a static library D2K3Pi
add_subdirectory("${CMAKE_SOURCE_DIR}/lib")

# Debug
set(CMAKE_BUILD_TYPE Debug)

# Bin ampgen decays
add_subdirectory("${CMAKE_SOURCE_DIR}/test/bin_decays/")

# UT
enable_testing()
add_subdirectory("${CMAKE_SOURCE_DIR}/test/ut/")

# IT
# current noop
add_subdirectory("${CMAKE_SOURCE_DIR}/test/it/")

# Stupid manual tests
add_subdirectory("${CMAKE_SOURCE_DIR}/test/manual/")

# Pull study
add_subdirectory("${CMAKE_SOURCE_DIR}/test/pull_study/")
