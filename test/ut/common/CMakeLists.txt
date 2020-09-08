set(commonUT
 ../driver.cpp
 test_util.cpp)

add_executable(ut_common.exe ${commonUT})

target_link_libraries(ut_common.exe common charmFitter ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_test(NAME commonUT COMMAND ut_common.exe)
