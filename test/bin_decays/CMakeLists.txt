add_executable(bin_decays.exe bin_generated_decays.cpp ReadAmpGen.cpp)
target_link_libraries(bin_decays.exe PRIVATE D2K3Pi PRIVATE ${ROOT_LIBRARIES} PRIVATE ${Boost_LIBRARIES} ${CMAKE_DL_LIBS})
