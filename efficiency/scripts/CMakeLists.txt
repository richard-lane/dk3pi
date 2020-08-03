# Dalitz plot thing
add_executable(dalitzthing.exe ${CMAKE_SOURCE_DIR}/efficiency/scripts/dalitzPlot.cpp)
target_link_libraries(dalitzthing.exe efficiency D2K3Pi ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
