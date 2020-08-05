# Dalitz plot thing
add_executable(dalitzthing.exe ${CMAKE_SOURCE_DIR}/efficiency/scripts/dalitzPlot.cpp)
target_link_libraries(dalitzthing.exe efficiency D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

# Weird test
add_executable(test_hist.exe ${CMAKE_SOURCE_DIR}/efficiency/scripts/test_plot_hist.cpp)
target_link_libraries(test_hist.exe efficiency D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
