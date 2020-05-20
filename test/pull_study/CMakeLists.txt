add_executable(pull-study pull_study.cpp PullStudyHelpers.cpp ${CMAKE_SOURCE_DIR}/test/pull_study/DecaySimulator.cpp)
target_link_libraries(pull-study D2K3Pi ${Boost_LIBRARIES} ${ROOT_LIBRARIES})
