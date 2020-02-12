# Add lib/ to our include path so that we can access headers from our UT
include_directories("${CMAKE_SOURCE_DIR}/lib/")
include_directories("${CMAKE_SOURCE_DIR}/src")

find_package(Boost REQUIRED COMPONENTS unit_test_framework filesystem system)

add_executable(manual_tests.exe test_simulator_equations.cpp)
target_link_libraries(manual_tests.exe D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_test(NAME manual COMMAND manual_tests.exe)
