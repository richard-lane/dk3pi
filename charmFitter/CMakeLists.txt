set(fitter
  src/fitter/CharmFitterBase.cpp
  src/fitter/PolynomialFitter.cpp
  src/fitter/ConstrainedFitter.cpp
  src/fitter/UnconstrainedFitter.cpp
)

add_library(charmFitter STATIC
  src/PhaseSpaceBinning.cpp
  src/RatioCalculator.cpp
  src/DecaySimulator.cpp
  src/physics.cpp
  src/fitterUtil.cpp
  ${fitter}
  )

target_include_directories(charmFitter
  PUBLIC "${CMAKE_SOURCE_DIR}/charmFitter/include"
  PUBLIC "${CMAKE_SOURCE_DIR}/charmFitter/include/fitter"
  PUBLIC "${CMAKE_SOURCE_DIR}/common/include"
  )

target_compile_features(charmFitter
  PUBLIC cxx_auto_type cxx_range_for cxx_uniform_initialization
  )

add_subdirectory(scripts)
