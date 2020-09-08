find_package(Boost REQUIRED COMPONENTS unit_test_framework filesystem system)

# Define main() fcn
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)

add_subdirectory(common/)
add_subdirectory(charmFitter/)
add_subdirectory(efficiency/)
