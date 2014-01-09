cmake_minimum_required(VERSION 2.6.2)

add_definitions(-std=c++0x -lconfig++)

set(BIOCPP_LIBRARY_FOLDER ${BIOCPP_LIBRARY_FOLDER} PARENT_SCOPE)

add_subdirectory(chemical_component_dictionary)
add_subdirectory(fasta)
add_subdirectory(proteins)
add_subdirectory(utils)
add_subdirectory(standard)

add_library(biocpp-core chemical_component_dictionary/atoms/atom_dictionary.cpp 
                        chemical_component_dictionary/elements/element_dictionary.cpp
                        chemical_component_dictionary/residues/residue_dictionary.cpp
                        fasta/NeedlemanWunsch.cpp 
                        fasta/StrictNeedlemanWunsch.cpp
                        proteins/dpss/dpss_id.cpp
                        proteins/dpss/h_bridge_energy.cpp
                        proteins/atom_position_reconstruction/bb_hydrogen.cpp
                        proteins/atom_position_reconstruction/beta_carbon.cpp
                        utils/errors_and_warnings/errors_and_warnings.cpp
           )
