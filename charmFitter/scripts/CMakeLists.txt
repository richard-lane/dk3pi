# Plots showing the behaviour of the decay simulator + fitters
add_executable(plots plots.cpp)
target_link_libraries(plots charmFitter common ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

# Pull study
add_executable(pull-study pull_study.cpp) 
target_link_libraries(pull-study common charmFitter ${Boost_LIBRARIES} ${ROOT_LIBRARIES})

# Test that our parameter scan works ok for a fit to a simple polynomial
add_executable(test_param_scan.exe test_param_scan.cpp)
target_link_libraries(test_param_scan.exe common charmFitter ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

# Example decay binning
add_executable(bin_decays.exe bin_generated_decays.cpp ReadAmpGen.cpp)
target_link_libraries(bin_decays.exe common charmFitter ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_DL_LIBS})

# Check the landscape of the CLEO likelihood
add_executable(cleo_likelihood.exe cleo_likelihood.cpp)
target_link_libraries(cleo_likelihood.exe common charmFitter ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

# Combine CLEO likelihood with toy data
add_executable(cleo_combination.exe cleo_combination.cpp)
target_link_libraries(cleo_combination.exe charmFitter common ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
