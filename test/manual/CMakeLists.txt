# Test that our parameter scan works ok for a fit to a simple polynomial
add_executable(test_param_scan.exe test_param_scan.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/PullStudyHelpers.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/DecaySimulator.cpp)
target_link_libraries(test_param_scan.exe common charmFitter ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
