find_package(Boost REQUIRED COMPONENTS unit_test_framework filesystem system)

add_executable(ut.exe ut.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/PullStudyHelpers.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/DecaySimulator.cpp)

target_link_libraries(ut.exe common charmFitter efficiency ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} dl)

#target_include_directories(ut.exe PUBLIC "${CMAKE_SOURCE_DIR}/common/include/" "${CMAKE_SOURCE_DIR}/charmFitter/include/" "${CMAKE_SOURCE_DIR}/src/" "${CMAKE_SOURCE_DIR}/efficiency/include")

add_test(NAME ut COMMAND ut.exe)
