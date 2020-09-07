add_executable(pull-study pull_study.cpp PullStudyHelpers.cpp DecaySimulator.cpp)

target_link_libraries(pull-study common charmFitter ${Boost_LIBRARIES} ${ROOT_LIBRARIES})

target_include_directories(pull-study PUBLIC "${CMAKE_SOURCE_DIR}/charmFitter/include/")
