set(fitterUT
  ../driver.cpp
  test_fitter.cpp
  test_fitter_utils.cpp
  test_pull_study.cpp
  test_ratio_calculator.cpp
  test_simulator.cpp)

add_executable(ut_fitter.exe
  ${fitterUT}
  ${CMAKE_SOURCE_DIR}/test/pull_study/PullStudyHelpers.cpp
  ${CMAKE_SOURCE_DIR}/test/pull_study/DecaySimulator.cpp)

target_link_libraries(ut_fitter.exe common charmFitter ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

target_include_directories(ut_fitter.exe PUBLIC "${CMAKE_SOURCE_DIR}/test/pull_study/")

add_test(NAME fitterUT COMMAND ut_fitter.exe)
