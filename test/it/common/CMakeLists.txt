set(commonIT
 ../driver.cpp
 test_correlated_random_numbers.cpp)

add_executable(it_common.exe ${commonIT})

target_link_libraries(it_common.exe common charmFitter ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_test(NAME commonIT COMMAND it_common.exe)
