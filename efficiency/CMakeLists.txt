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

# Weird AmpGen amplitude model stuff
set(ampGen
  ampGen/src/findZ.cpp
  ampGen/src/flatPhsp.cpp
  ampGen/src/amplitudes.cpp)

# Define monolithic efficiency library
add_library(efficiency STATIC
  ${misc}
  ${chowLiu}
  ${ampGen}
  )

# Include dirs
set(efficiencyIncludeDirs
  "misc/include/"
  "chowLiu/include/" 
  "ampGen/include/"
  "${CMAKE_SOURCE_DIR}/common/include")
target_include_directories(efficiency
  PUBLIC ${efficiencyIncludeDirs}
  )

target_compile_features(efficiency
  PUBLIC cxx_auto_type cxx_range_for cxx_uniform_initialization
  )

add_subdirectory(scripts/)
