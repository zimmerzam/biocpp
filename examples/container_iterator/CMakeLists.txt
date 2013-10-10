
project( container_iterator )
set_property(
  SOURCE ${PROJECT_NAME}.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_PDB 
)
add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)
