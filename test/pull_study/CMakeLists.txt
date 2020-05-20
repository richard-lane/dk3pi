add_executable(pull-study pull_study.cpp PullStudyHelpers.cpp)
target_link_libraries(pull-study D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
