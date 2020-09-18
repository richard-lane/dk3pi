# I guess we need Python for this library
# 3.7 is the version I'm using, so use that
find_package(Python3 3.7 REQUIRED COMPONENTS Interpreter Development NumPy)
message(STATUS "Python3 executable: " ${Python3_EXECUTABLE})

if (NOT DEFINED Python3_INCLUDE_DIRS OR NOT DEFINED Python3_NumPy_INCLUDE_DIRS)
  message(SEND_ERROR "Python include dirs not set")
else()
  message(STATUS "Python3 include dir: " ${Python3_INCLUDE_DIRS})
endif(NOT DEFINED Python3_INCLUDE_DIRS OR NOT DEFINED Python3_NumPy_INCLUDE_DIRS)

# Miscellaneous utility stuff
set(misc
  misc/src/efficiencyUtil.cpp
  misc/src/ReadRoot.cpp)
  
# Chow Liu efficiency estimate
set(chowLiu
  chowLiu/src/graphTheory.cpp
  chowLiu/src/ChowLiu.cpp
  chowLiu/src/efficiency.cpp
  )

# BDT efficiency estimate
set(bdt
  bdt_reweighting/src/bdt_reweighting.cpp
  )

# Weird AmpGen amplitude model stuff
set(ampGen
  ampGen/src/findZ.cpp
  ampGen/src/flatPhsp.cpp
  ampGen/src/amplitudes.cpp)

# Define monolithic efficiency library
add_library(efficiency STATIC
  ${misc}
  ${chowLiu}
  ${bdt}
  ${ampGen}
  )

target_link_libraries(efficiency PUBLIC ${Python3_LIBRARIES})

# Include dirs
set(efficiencyIncludeDirs
  "misc/include/"
  "chowLiu/include/" 
  "bdt_reweighting/include/"
  "ampGen/include/"
  "${CMAKE_SOURCE_DIR}/common/include"
  "${Python3_INCLUDE_DIRS}"
  "${Python3_NumPy_INCLUDE_DIRS}")
target_include_directories(efficiency
  PUBLIC ${efficiencyIncludeDirs}
  )

# Python path
# This only affects the path during CMake's execution, but we can make CTest tests inherit it
set(ENV{PYTHONPATH} $ENV{PYTHONPATH}:${CMAKE_CURRENT_SOURCE_DIR}/bdt_reweighting/)

target_compile_features(efficiency
  PUBLIC cxx_auto_type cxx_range_for cxx_uniform_initialization
  )

add_subdirectory(scripts/)
