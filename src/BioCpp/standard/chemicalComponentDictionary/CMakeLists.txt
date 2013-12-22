set(BIOCPP_LIBRARY_FOLDER ${BIOCPP_LIBRARY_FOLDER} PARENT_SCOPE)
set(LIBRARY_OUTPUT_PATH ${BIOCPP_LIBRARY_FOLDER})

add_library(biocpp-std-ccd element_dictionary_standard.cpp
                           atom_dictionary_standard.cpp
                           residue_dictionary_standard.cpp
           )
