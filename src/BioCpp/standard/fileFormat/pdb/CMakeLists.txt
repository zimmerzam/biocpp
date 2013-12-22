set(BIOCPP_LIBRARY_FOLDER ${BIOCPP_LIBRARY_FOLDER} PARENT_SCOPE)
set(LIBRARY_OUTPUT_PATH ${BIOCPP_LIBRARY_FOLDER})

add_library(biocpp-std-io sections_and_records/sections_and_records.cpp
                          ../../chemicalComponentDictionary/element_dictionary_standard.cpp
                          ../../chemicalComponentDictionary/atom_dictionary_standard.cpp
                          ../../chemicalComponentDictionary/residue_dictionary_standard.cpp
                          file_reader.cpp 
           )
