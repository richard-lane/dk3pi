# Add lib/ to our include path so that we can access headers from our UT
include_directories("${CMAKE_SOURCE_DIR}/lib/")
include_directories("${CMAKE_SOURCE_DIR}/src")

find_package(Boost REQUIRED COMPONENTS unit_test_framework filesystem system)

add_executable(ut.exe ut.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/PullStudyHelpers.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/DecaySimulator.cpp)
target_link_libraries(ut.exe D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_test(NAME ut COMMAND ut.exe)
