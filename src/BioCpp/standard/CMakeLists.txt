set(BIOCPP_LIBRARY_FOLDER ${BIOCPP_LIBRARY_FOLDER} PARENT_SCOPE)
set(LIBRARY_OUTPUT_PATH ${BIOCPP_LIBRARY_FOLDER})

add_subdirectory(alignmentMatrix)
add_subdirectory(chemicalComponentDictionary)
add_subdirectory(fileFormat)
add_subdirectory(proteinDefinitions)

add_library(biocpp-std     
                           alignmentMatrix/BLOSUM45.cpp
                           alignmentMatrix/BLOSUM62.cpp
                           alignmentMatrix/BLOSUM80.cpp
                           alignmentMatrix/PAM30.cpp
                           alignmentMatrix/PAM70.cpp
                           alignmentMatrix/ZIMM1.cpp
                           
                           ../chemical_component_dictionary/atoms/atom_dictionary.cpp 
                           ../chemical_component_dictionary/elements/element_dictionary.cpp
                           ../chemical_component_dictionary/residues/residue_dictionary.cpp
                           chemicalComponentDictionary/element_dictionary_standard.cpp
                           chemicalComponentDictionary/atom_dictionary_standard.cpp
                           chemicalComponentDictionary/residue_dictionary_standard.cpp
                           
                           fileFormat/pdb/sections_and_records/sections_and_records.cpp
                           fileFormat/pdb/file_reader.cpp
                           
                           proteinDefinitions/base/dpss/h_bridge_map_constructor.cpp
           )
