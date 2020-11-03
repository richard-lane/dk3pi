# Dalitz plot thing
add_executable(dalitzthing.exe ${CMAKE_SOURCE_DIR}/efficiency/scripts/dalitzPlot.cpp)
target_link_libraries(dalitzthing.exe common efficiency charmFitter ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

# Example toy efficiency thing with 1d histogram projections too
add_executable(toy_efficiency_example.exe ${CMAKE_CURRENT_SOURCE_DIR}/scriptUtils.cpp ${CMAKE_CURRENT_SOURCE_DIR}/toy_efficiency_example.cpp)
target_link_libraries(toy_efficiency_example.exe common efficiency charmFitter ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

# Example reweighting prompt to semileptonic
add_executable(prompt_to_semileptonic.exe ${CMAKE_CURRENT_SOURCE_DIR}/prompt_to_semileptonic.cpp scriptUtils.cpp)
target_link_libraries(prompt_to_semileptonic.exe common efficiency ${ROOT_LIBRARIES})
