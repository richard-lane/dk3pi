# Pull study
add_executable(pull-study pull_study.cpp) 
target_link_libraries(pull-study common charmFitter ${Boost_LIBRARIES} ${ROOT_LIBRARIES})

# Test that our parameter scan works ok for a fit to a simple polynomial
add_executable(test_param_scan.exe
  test_param_scan.cpp)
target_link_libraries(test_param_scan.exe common charmFitter ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

# Example decay binning
add_executable(bin_decays.exe bin_generated_decays.cpp ReadAmpGen.cpp)
target_link_libraries(bin_decays.exe common charmFitter ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_DL_LIBS})
