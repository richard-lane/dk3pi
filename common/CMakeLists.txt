add_library(common STATIC
  ${CMAKE_SOURCE_DIR}/common/lib/util.cpp
  )

target_include_directories(common PUBLIC include/)
