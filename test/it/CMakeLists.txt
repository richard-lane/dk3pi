find_package(Boost REQUIRED COMPONENTS unit_test_framework filesystem system)

add_executable(it.exe
               it.cpp
               test_correlated_random_numbers.cpp
               test_minuit_fitter.cpp
               test_ratio_calculator.cpp
               test_simulator.cpp
               efficiency/test_amplitudes.cpp
               # efficiency/flat_phase_space.cpp  # What the heck
               efficiency/test_chow_liu.cpp
               ${CMAKE_SOURCE_DIR}/test/pull_study/DecaySimulator.cpp)

target_include_directories(it.exe PUBLIC "${CMAKE_SOURCE_DIR}/charmFitter/include/" "${CMAKE_SOURCE_DIR}/efficiency/include")

target_link_libraries(it.exe common charmFitter efficiency ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} dl)

add_test(NAME it COMMAND it.exe)
