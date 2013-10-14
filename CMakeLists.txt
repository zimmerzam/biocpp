cmake_minimum_required(VERSION 2.6.2)

## Policy settings ##
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

## store these variables cause CGAL usually overrides them ##
set( BUILD_TYPE ${CMAKE_BUILD_TYPE} )
set( CXX_FLAGS "${CMAKE_CXX_FLAGS}" )

## Find Boost and its components ##
find_package( Boost REQUIRED )
if( Boost_FOUND )
  include_directories( ${Boost_INCLUDE_DIRS} )
else (  )
  message(STATUS "Projects that require the Boost library will not be compiled.")
endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )
if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
else (  )
  message(STATUS "Projects that require the CGAL library will not be compiled.")
endif()

if( CXX=icc OR CXX=icpc )
  set( CMAKE_CXX_FLAGS "-fp-model strict" )
endif()

# Set compiler flags and build type
set( CMAKE_BUILD_TYPE ${BUILD_TYPE} )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_FLAGS} -std=c++0x " )
set( CMAKE_CXX_FLAGS_DEBUG " -Wall")
set( CMAKE_CXX_FLAGS_RELEASE " -msse2 -O2")

## This is the project base folder ##
set( BIOCPP_BASE ${CMAKE_SOURCE_DIR} )

## Configure version number ##
set( BIOCPP_VERSION_MAJOR 1 )
set( BIOCPP_VERSION_MINOR 0 )
configure_file(
  "${BIOCPP_BASE}/src/BioCpp/version.h.ini"
  "${BIOCPP_BASE}/src/BioCpp/version.h"
)

## Add rules to build all the examples ##
set( BIOCPP_EXAMPLES ${BIOCPP_BASE}/examples )
add_subdirectory( ${BIOCPP_EXAMPLES} )
