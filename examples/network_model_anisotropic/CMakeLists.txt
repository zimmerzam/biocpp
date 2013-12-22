
project( network_model_anisotropic )
set_property(
  SOURCE ${PROJECT_NAME}.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_PDB
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_STANDARD
)
add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)
target_link_libraries(${PROJECT_NAME} pdb)
