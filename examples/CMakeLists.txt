cmake_minimum_required(VERSION 2.6.2)

add_definitions( -msse2 -O2 -Wall -std=c++0x -static )

project( container )
add_executable(container container.cpp)

project( ids )
set_property(
  SOURCE ids.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_ID 
)
add_executable(ids ids.cpp)

project( pdb )
set_property(
  SOURCE pdb.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_PDB 
)
add_executable(pdb pdb.cpp)

project( fasta_align )
set_property(
  SOURCE fasta_align.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_FASTA_ALIGNMENT
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_BLOSUM62 
)
add_executable(fasta_align fasta_align.cpp)

project( strict_align )
set_property(
  SOURCE strict_align.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_FASTA_ALIGNMENT
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_FASTA_SUBST_MATRIX_ZIMM1 
)
add_executable(strict_align strict_align.cpp)

project( iterate )
set_property(
  SOURCE iterate.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_PDB 
)
add_executable(iterate iterate.cpp)

project( dpss )
set_property(
  SOURCE dpss.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_PDB
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_DPSS
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_STANDARD 
)
add_executable(dpss dpss.cpp)

project( surface_area_lcpo )
set_property(
  SOURCE surface_area_lcpo.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_PDB
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_MORPHOLOGY
  PROPERTY COMPILE_DEFINITIONS BIOCPP_DONT_INCLUDE_CGAL
)
add_executable(surface_area_lcpo surface_area_lcpo.cpp)

project( contact_map )
set_property(
  SOURCE contact_map.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_STANDARD
)
add_executable(contact_map contact_map.cpp)

project( gaussian_network_model )
add_executable(gaussian_network_model gaussian_network_model.cpp)

project( anisotropic_network_model )
set_property(
  SOURCE anisotropic_network_model.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_STANDARD
)
add_executable(anisotropic_network_model anisotropic_network_model.cpp)

project( mutations )
set_property(
  SOURCE mutations.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_STANDARD
)
add_executable(mutations mutations.cpp)

project( C_alpha_variance )
set_property(
  SOURCE C_alpha_variance.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_STANDARD
)
add_executable(C_alpha_variance C_alpha_variance.cpp)

project( print_atom_coordinate )
set_property(
  SOURCE print_atom_coordinate.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_STANDARD
)
add_executable(print_atom_coordinate print_atom_coordinate.cpp)

project( complex_constructor )
set_property(
  SOURCE complex_constructor.cpp
  PROPERTY COMPILE_DEFINITIONS BIOCPP_INCLUDE_STANDARD
)
add_executable(complex_constructor complex_constructor.cpp)

project( has_coordinate )
add_executable(has_coordinate has_coordinate.cpp)
