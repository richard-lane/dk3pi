add_executable(example_s_weighting.exe ${CMAKE_SOURCE_DIR}/lhcbMonteCarlo/sWeighting/scripts/example_s_weighting.cpp)
target_link_libraries(example_s_weighting.exe sWeighting ${ROOT_LIBRARIES})
