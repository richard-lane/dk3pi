add_executable(add_delta_m_branch.exe ${CMAKE_SOURCE_DIR}/lhcbMonteCarlo/sWeighting/scripts/add_delta_m_branch.cpp)
target_include_directories(add_delta_m_branch.exe PUBLIC ${CMAKE_SOURCE_DIR}/common/include/)
target_link_libraries(add_delta_m_branch.exe ${ROOT_LIBRARIES} stdc++fs)

add_executable(example_s_weighting.exe ${CMAKE_SOURCE_DIR}/lhcbMonteCarlo/sWeighting/scripts/example_s_weighting.cpp)
target_link_libraries(example_s_weighting.exe sWeighting ${ROOT_LIBRARIES} common)

add_executable(make_cuts.exe ${CMAKE_SOURCE_DIR}/lhcbMonteCarlo/sWeighting/scripts/makeCuts.cpp)
target_link_libraries(make_cuts.exe cuts ${ROOT_LIBRARIES})
