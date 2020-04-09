set(fitter
  fitter/Fitter.cpp
  fitter/FitterUtils.cpp
  fitter/BaseFitter.cpp
  fitter/RootFitter.cpp
  fitter/MinuitFitterBase.cpp
  fitter/MinuitPolynomialFitter.cpp
  fitter/MinuitPolynomialScan.cpp
  )

add_library(D2K3Pi SHARED
  DecaySimulator.cpp
  ${fitter}
  MinuitFitter.cpp
  PhaseSpaceBinning.cpp
  RatioCalculator.cpp
  ReadAmpGen.cpp
  util.cpp
  )

set(CMAKE_BUILD_TYPE Debug)

target_include_directories(D2K3Pi
  PUBLIC ${CMAKE_CURRENT_LIST_DIR} "../include/" "../src/" "fitter/"
  )

target_compile_features(D2K3Pi
  PUBLIC cxx_auto_type cxx_range_for cxx_uniform_initialization
  )