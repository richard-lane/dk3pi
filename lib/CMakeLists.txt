set(fitter
  src/fitter/FitterUtils.cpp
  src/fitter/BaseFitter.cpp
  src/fitter/RootFitter.cpp
  src/fitter/MinuitFitterBase.cpp
  src/fitter/MinuitScannerBase.cpp
  src/fitter/MinuitPolynomialFitter.cpp
  src/fitter/PhysicalFitter.cpp
  )

add_library(D2K3Pi STATIC
  src/PhaseSpaceBinning.cpp
  src/RatioCalculator.cpp
  src/util.cpp
  ${fitter}
  )

set(CMAKE_BUILD_TYPE Debug)

target_include_directories(D2K3Pi
  PUBLIC ${CMAKE_CURRENT_LIST_DIR} "include/" "include/fitter/"
  )

target_compile_features(D2K3Pi
  PUBLIC cxx_auto_type cxx_range_for cxx_uniform_initialization
  )
