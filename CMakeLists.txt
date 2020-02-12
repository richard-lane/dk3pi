cmake_minimum_required(VERSION 3.2)
project(D2K3Pi VERSION 0.1.0)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ROOT enforces C++14
set(CMAKE_CXX_STANDARD 14)

# ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include(${ROOT_USE_FILE})

# Boost
find_package( Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Add our core library
# This defines a static library D2K3Pi
add_subdirectory("${CMAKE_SOURCE_DIR}/lib")

# Some standard flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wformat-security -Werror -O0 -g")

# Pull study
add_executable(pull.exe src/pull_study.cpp src/PullStudyHelpers.cpp)
target_link_libraries(pull.exe PRIVATE D2K3Pi PRIVATE ${ROOT_LIBRARIES} PRIVATE ${Boost_LIBRARIES})

# Bin AmpGen's generated decays
add_executable(bin_decays.exe src/bin_generated_decays.cpp)
target_link_libraries(bin_decays.exe PRIVATE D2K3Pi PRIVATE ${ROOT_LIBRARIES} PRIVATE ${Boost_LIBRARIES} ${CMAKE_DL_LIBS})

# UT
enable_testing()
add_subdirectory("${CMAKE_SOURCE_DIR}/test/ut/")

# IT
# current noop
add_subdirectory("${CMAKE_SOURCE_DIR}/test/it/")

# Stupid manual tests
add_subdirectory("${CMAKE_SOURCE_DIR}/test/manual/")
