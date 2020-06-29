add_library(efficiency STATIC
  lib/amplitudes.cpp
  )

set(CMAKE_BUILD_TYPE Debug)

target_include_directories(efficiency
  PUBLIC ${CMAKE_CURRENT_LIST_DIR} "include/" "include/"
  )

target_compile_features(efficiency
  PUBLIC cxx_auto_type cxx_range_for cxx_uniform_initialization
  )
