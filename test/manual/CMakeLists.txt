# Add lib/ to our include path so that we can access headers from our UT
include_directories("${CMAKE_SOURCE_DIR}/lib/")
include_directories("${CMAKE_SOURCE_DIR}/src/")

find_package(Boost REQUIRED COMPONENTS unit_test_framework filesystem system)

# Simulator equations: check that we're generating sensible looking decay times with our accept-reject
add_executable(simulator_equations.exe ${CMAKE_SOURCE_DIR}/src/PullStudyHelpers.cpp test_simulator_equations.cpp)
target_link_libraries(simulator_equations.exe D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
add_test(NAME simulator_equations COMMAND simulator_equations.exe)

# Fitter: test that our fitter doesn't do anything weird
add_executable(test_fitter.exe test_fitter.cpp ${CMAKE_SOURCE_DIR}/src/PullStudyHelpers.cpp)
target_include_directories(test_fitter.exe PUBLIC ${CMAKE_SOURCE_DIR}/src/})
target_link_libraries(test_fitter.exe D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
add_test(NAME test_fitter COMMAND test_fitter.exe)
