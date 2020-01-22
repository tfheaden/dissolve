add_library(ff
  # Canongia-Lopes, Deschamps, & Padua IL 2010
  cldp-il-2010.cpp
  cldp-il-2010.h

  # SPC/Fw
  spcfw.cpp
  spcfw.h

  # UFF
  uff.cpp
  uffatomtype.cpp
  uff.h
  uffatomtype.h
)

add_subdirectory(ludwig)
add_subdirectory(oplsaa2005)

include_directories(
../
${CMAKE_SOURCE_DIR}
${CMAKE_BINARY_DIR}
${CMAKE_BINARY_DIR}/src
)
