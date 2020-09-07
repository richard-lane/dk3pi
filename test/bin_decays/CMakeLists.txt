add_executable(bin_decays.exe bin_generated_decays.cpp ReadAmpGen.cpp)
target_link_libraries(bin_decays.exe common charmFitter ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_DL_LIBS})
