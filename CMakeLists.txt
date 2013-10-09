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

# store these variables cause CGAL will override them..
SET( BUILD_TYPE ${CMAKE_BUILD_TYPE} )
SET( CXX_FLAGS "${CMAKE_CXX_FLAGS}" )

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )
if ( CGAL_FOUND )
  # include helper file
  include( ${CGAL_USE_FILE} )
endif()
if ( NOT CGAL_FOUND )
  message(STATUS "Projects that require the CGAL library will not be compiled.")
  return()  
endif()

# Boost and its components
find_package( Boost REQUIRED )
if ( NOT Boost_FOUND )
  message(STATUS "Projects that require the Boost library will not be compiled.")
endif()

# restore these variables to its original value
SET( CMAKE_BUILD_TYPE ${BUILD_TYPE} )
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_FLAGS}" )

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x") 
SET( CMAKE_CXX_FLAGS_DEBUG " -Wall")
SET( CMAKE_CXX_FLAGS_RELEASE " -msse2 -O2")

SET( BIOCPP_BASE ${CMAKE_SOURCE_DIR} )
SET( BIOCPP_SOURCES ${BIOCPP_BASE}/src/BioCpp/ )
SET( BIOCPP_CONFIG_FOLDER ${BIOCPP_BASE}/cmake/BioCpp )

SET( BIOCPP_EXAMPLE_SOURCES ${BIOCPP_BASE}/examples )
SET( BIOCPP_EXAMPLE_CONFIG_FOLDER ${BIOCPP_BASE}/cmake/examples )

SET( BIOCPP_VERSION_MAJOR 1 )
SET( BIOCPP_VERSION_MINOR 0 )
CONFIGURE_FILE(
  "${BIOCPP_CONFIG_FOLDER}/version.h.ini"
  "${BIOCPP_BASE}/src/BioCpp/version.h"
)

SET( CMAKE_MODULE_PATH ${BIOCPP_BASE}/cmake/modules )

INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/lcpo_exposed.cmake )
INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/surface_area_lcpo.cmake )
INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/container.cmake )
INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/dpss.cmake )
INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/ids.cmake )
INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/iterate.cmake )
INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/pdb.cmake )
INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/fasta_align.cmake )
INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/strict_align.cmake )
INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/contact_map.cmake )
INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/gaussian_network_model.cmake )
INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/anisotropic_network_model.cmake )
INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/mutations.cmake )
INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/C_alpha_variance.cmake )
INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/print_atom_coordinate.cmake )
INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/complex_constructor.cmake )
INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/has_coordinate.cmake )
INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/md_gas_of_identical_particles.cmake )

if( CGAL_FOUND AND Boost_FOUND )
  INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/alpha_exposed.cmake )
  INCLUDE( ${BIOCPP_EXAMPLE_CONFIG_FOLDER}/delaunay_exposed.cmake )
endif()
