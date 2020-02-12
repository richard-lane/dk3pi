add_library(D2K3Pi STATIC
  DecaySimulator.cpp
  Fitter.cpp
  MCGenerator.cpp
  PhaseSpaceBinning.cpp
  RatioCalculator.cpp
  ReadAmpGen.cpp
  util.cpp
  )

target_include_directories(D2K3Pi
  PUBLIC ${CMAKE_CURRENT_LIST_DIR} "../include/"
  )
target_compile_features(D2K3Pi
  PUBLIC cxx_auto_type cxx_range_for cxx_uniform_initialization
  )