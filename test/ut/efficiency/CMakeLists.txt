set(efficiencyUT
 ../driver.cpp
 test_amplitude.cpp
 test_efficiency.cpp
 test_efficiency_util.cpp
 test_flat_phsp.cpp
 test_graph_theory.cpp
 test_readroot.cpp
 test_z.cpp)

add_executable(ut_efficiency.exe ${efficiencyUT})

target_link_libraries(ut_efficiency.exe common efficiency ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} dl)

add_test(NAME efficiencyUT COMMAND ut_efficiency.exe)
