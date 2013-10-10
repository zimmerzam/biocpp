
project( ids )
set_property(
  SOURCE ${PROJECT_NAME}.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_ID 
)
add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)
