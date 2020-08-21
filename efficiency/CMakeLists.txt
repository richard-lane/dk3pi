add_library(efficiency STATIC
  ${CMAKE_SOURCE_DIR}/lib/src/util.cpp
  lib/amplitudes.cpp
  lib/ReadRoot.cpp
  lib/findZ.cpp
  lib/flatPhsp.cpp
  lib/efficiency.cpp
  lib/graphTheory.cpp
  lib/efficiencyUtil.cpp
  )

set(CMAKE_BUILD_TYPE Debug)

target_include_directories(efficiency
  PUBLIC ${CMAKE_CURRENT_LIST_DIR} "include/" "include/" "${CMAKE_SOURCE_DIR}/lib/include/"
  )

target_compile_features(efficiency
  PUBLIC cxx_auto_type cxx_range_for cxx_uniform_initialization
  )
