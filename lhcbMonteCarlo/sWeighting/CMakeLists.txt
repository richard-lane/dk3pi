# sWeighting lib
add_library(sWeighting STATIC
  ${CMAKE_SOURCE_DIR}/lhcbMonteCarlo/sWeighting/sWeighting.cpp
  )
target_link_libraries(sWeighting stdc++fs) # Filesystem lib

include_directories(${CMAKE_SOURCE_DIR}/lhcbMonteCarlo/sWeighting/)
add_subdirectory("${CMAKE_SOURCE_DIR}/lhcbMonteCarlo/sWeighting/scripts/")

# cuts lib
add_library(cuts STATIC ${CMAKE_SOURCE_DIR}/lhcbMonteCarlo/sWeighting/cuts.cpp)
