add_executable(example_s_weighting.exe ${CMAKE_SOURCE_DIR}/lhcbMonteCarlo/sWeighting/scripts/example_s_weighting.cpp)
target_include_directories(example_s_weighting.exe PUBLIC ${CMAKE_SOURCE_DIR/common/include/})
target_link_libraries(example_s_weighting.exe sWeighting ${ROOT_LIBRARIES} common)
