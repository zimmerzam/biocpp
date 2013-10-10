
project( contact_map )
set_property(
  SOURCE ${PROJECT_NAME}.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_STANDARD
)
add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)
