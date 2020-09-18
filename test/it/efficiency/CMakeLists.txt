set(efficiencyIT
 ../driver.cpp
 flat_phase_space.cpp
 test_amplitudes.cpp
 test_chow_liu.cpp
 test_bdt.cpp)

add_executable(it_efficiency.exe ${efficiencyIT})

target_link_libraries(it_efficiency.exe PUBLIC common efficiency ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} dl)

add_test(NAME efficiencyIT COMMAND it_efficiency.exe)

add_test(NAME efficiencyPythonBDT COMMAND ${Python3_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR} -sv)
set_tests_properties(efficiencyPythonBDT PROPERTIES ENVIRONMENT PYTHONPATH=$ENV{PYTHONPATH})
