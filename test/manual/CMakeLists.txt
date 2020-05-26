# Add lib/ to our include path so that we can access headers from our UT
include_directories("${CMAKE_SOURCE_DIR}/lib/")
include_directories("${CMAKE_SOURCE_DIR}/src/")
include_directories("${CMAKE_SOURCE_DIR}/lib/fitter/")

find_package(Boost REQUIRED COMPONENTS unit_test_framework filesystem system)

# Simulator equations: check that we're generating sensible looking decay times with our accept-reject
add_executable(simulator_equations.exe ${CMAKE_SOURCE_DIR}/test/pull_study/PullStudyHelpers.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/DecaySimulator.cpp test_simulator_equations.cpp)
target_link_libraries(simulator_equations.exe D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

# Test that our minuit fitter and the ROOT builtin fitter do sensible things
add_executable(test_minuit_fitter.exe test_minuit_fitters.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/PullStudyHelpers.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/DecaySimulator.cpp)
target_link_libraries(test_minuit_fitter.exe D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

# Test that our parameter scan works ok for a fit to a simple polynomial
add_executable(test_param_scan.exe test_param_scan.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/PullStudyHelpers.cpp testutil.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/DecaySimulator.cpp)
target_link_libraries(test_param_scan.exe D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

# Pull study on generated decay times, fitting to  (a + bt + ct^2)exp(-width * t)
add_executable(test_simulator_pull.exe SimulatorEqnsPull.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/PullStudyHelpers.cpp testutil.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/DecaySimulator.cpp)
target_link_libraries(test_simulator_pull.exe D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

# Pull study on ideal ratios with gaussian noise
add_executable(test_fitter.exe pull_fitter.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/PullStudyHelpers.cpp testutil.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/DecaySimulator.cpp)
target_link_libraries(test_fitter.exe D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

# Pull study generating counts by integrating the expected rates
add_executable(test_fitter_calculator.exe pull_calculator_and_fitter.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/PullStudyHelpers.cpp testutil.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/DecaySimulator.cpp)
target_link_libraries(test_fitter_calculator.exe D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
