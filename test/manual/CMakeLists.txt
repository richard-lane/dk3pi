# Add lib/ to our include path so that we can access headers from our UT
include_directories("${CMAKE_SOURCE_DIR}/lib/")
include_directories("${CMAKE_SOURCE_DIR}/src/")
include_directories("${CMAKE_SOURCE_DIR}/lib/fitter/")

find_package(Boost REQUIRED COMPONENTS unit_test_framework filesystem system)

# Simulator equations: check that we're generating sensible looking decay times with our accept-reject
add_executable(simulator_equations.exe ${CMAKE_SOURCE_DIR}/src/PullStudyHelpers.cpp test_simulator_equations.cpp)
target_link_libraries(simulator_equations.exe D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
add_test(NAME simulator_equations COMMAND simulator_equations.exe)

# Fitter: test that our fitter doesn't do anything weird
add_executable(test_fitter.exe test_fitter.cpp ${CMAKE_SOURCE_DIR}/src/PullStudyHelpers.cpp testutil.cpp)
target_link_libraries(test_fitter.exe D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
add_test(NAME test_fitter COMMAND test_fitter.exe)

# Test that our minuit fitter and the ROOT builtin fitter do sensible things
add_executable(test_minuit_fitter.exe test_minuit_fitters.cpp ${CMAKE_SOURCE_DIR}/src/PullStudyHelpers.cpp)
target_link_libraries(test_minuit_fitter.exe D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
add_test(NAME test_compare_minuit_root COMMAND test_minuit_fitter.exe)

# Test that our parameter scan works ok for a fit to a simple polynomial
add_executable(test_param_scan.exe test_param_scan.cpp ${CMAKE_SOURCE_DIR}/src/PullStudyHelpers.cpp testutil.cpp)
target_link_libraries(test_param_scan.exe D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
add_test(NAME test_param_scan COMMAND test_param_scan.exe)

# Pull study
# This is so slow that I won't add a make-test to it
add_executable(pull-study test_pull_study.cpp ${CMAKE_SOURCE_DIR}/src/PullStudyHelpers.cpp)
target_link_libraries(pull-study D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
