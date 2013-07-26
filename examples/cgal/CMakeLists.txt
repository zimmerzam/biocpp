cmake_minimum_required(VERSION 2.6.2)

if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )
if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()

# include helper file
include( ${CGAL_USE_FILE} )

# Boost and its components
find_package( Boost REQUIRED )
if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()  
endif()

add_definitions( -msse2 -O2 -Wall -std=c++0x -static )

# Alpha Shape
project( alpha_shapes )
add_executable(alpha_shapes alpha_shapes.cpp)
target_link_libraries(alpha_shapes ${CGAL_LIBRARIES})

# Alpha Shape
project( triang3_boost_zip )
add_executable(triang3_boost_zip triang3_boost_zip.cpp)
target_link_libraries(triang3_boost_zip ${CGAL_LIBRARIES})

project( protein_surface )
set_property(
  SOURCE protein_surface.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_PDB
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_STANDARD
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_CGAL
)
add_executable(protein_surface protein_surface.cpp)
target_link_libraries(protein_surface ${CGAL_LIBRARIES})

# kernel
project( kernel )
set_property(
  SOURCE kernel.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_PDB 
)
add_executable(kernel kernel.cpp)
target_link_libraries(kernel ${CGAL_LIBRARIES})
