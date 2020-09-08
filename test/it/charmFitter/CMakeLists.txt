set(fitterIT
 ../driver.cpp
 test_minuit_fitter.cpp
 test_ratio_calculator.cpp
 test_simulator.cpp
 ${CMAKE_SOURCE_DIR}/test/pull_study/DecaySimulator.cpp)

add_executable(it_fitter.exe ${fitterIT})

target_link_libraries(it_fitter.exe common charmFitter ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_include_directories(it_fitter.exe PUBLIC "${CMAKE_SOURCE_DIR}/test/pull_study/")

add_test(NAME fitterIT COMMAND it_fitter.exe)
