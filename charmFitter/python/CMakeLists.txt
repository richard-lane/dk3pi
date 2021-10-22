# Library for a nice python-compatible interface with the CLEO scan
add_library(pyFitter SHARED
    ${CMAKE_SOURCE_DIR}/charmFitter/python/src/cleoFitter.cpp
    ${CMAKE_SOURCE_DIR}/charmFitter/python/src/simulator.cpp
    ${CMAKE_SOURCE_DIR}/charmFitter/python/src/bindings.cpp
)

find_package(Python COMPONENTS Interpreter Development)
find_package(PythonLibs REQUIRED)

# Embarrassing hack: find the include dir containing pybind11 with the below shell command
# I think the LCG envs are making python confusing and not automatically finding the right include dir
# Need to strip whitespace from the output in order to get sensible makefiles
execute_process(COMMAND bash -c "${Python_EXECUTABLE} -m pybind11 --includes"
                OUTPUT_VARIABLE PYBIND_INCLUDES
                OUTPUT_STRIP_TRAILING_WHITESPACE)

# Also need to remove "quotes" from the output; do this here
separate_arguments(PYBIND_INCLUDES UNIX_COMMAND "${PYBIND_INCLUDES}")

target_include_directories(pyFitter
  PUBLIC "${CMAKE_SOURCE_DIR}/charmFitter/CLEO"
  PUBLIC "${CMAKE_SOURCE_DIR}/charmFitter/include"
  PUBLIC "${CMAKE_SOURCE_DIR}/charmFitter/include/fitter"
  PUBLIC "${CMAKE_SOURCE_DIR}/charmFitter/python/include/"
  PUBLIC "${CMAKE_SOURCE_DIR}/common/include"
  PUBLIC "${PYTHON_INCLUDE_DIRS}"
  )

target_link_libraries(pyFitter charmFitter common ${PYTHON_LIBRARIES} ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

target_compile_features(pyFitter
  PUBLIC cxx_auto_type cxx_range_for cxx_uniform_initialization
  )
target_compile_options(pyFitter PRIVATE "${PYBIND_INCLUDES}")
execute_process(COMMAND bash -c "python3-config --extension-suffix"
                OUTPUT_VARIABLE LIBSUFFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE)
set_target_properties(pyFitter PROPERTIES OUTPUT_NAME "cleoScan"
                                          SUFFIX "${LIBSUFFIX}"
                                          PREFIX "lib")

# The python binding

