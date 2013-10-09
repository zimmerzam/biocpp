SET( BIOCPP_EXAMPLE_SOURCES ${BIOCPP_EXAMPLE_SOURCES} PARENT_SCOPE )

project( dpss )
set_property(
  SOURCE ${BIOCPP_EXAMPLE_SOURCES}/${PROJECT_NAME}.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_PDB
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_DPSS
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_STANDARD 
)
add_executable(${PROJECT_NAME} ${BIOCPP_EXAMPLE_SOURCES}/${PROJECT_NAME}.cpp)