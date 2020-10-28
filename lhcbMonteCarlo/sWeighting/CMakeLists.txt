add_library(sWeighting STATIC
  ${CMAKE_SOURCE_DIR}/lhcbMonteCarlo/sWeighting/sWeighting.cpp
  )

# Link filesystem lib
target_link_libraries(sWeighting stdc++fs)

include_directories(${CMAKE_SOURCE_DIR}/lhcbMonteCarlo/sWeighting/)

add_subdirectory("${CMAKE_SOURCE_DIR}/lhcbMonteCarlo/sWeighting/scripts/")
