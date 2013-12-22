set(BIOCPP_LIBRARY_FOLDER ${BIOCPP_LIBRARY_FOLDER} PARENT_SCOPE)
set(LIBRARY_OUTPUT_PATH ${BIOCPP_LIBRARY_FOLDER})

add_library(biocpp-ccd atoms/atom_dictionary.cpp 
                       elements/element_dictionary.cpp
                       residues/residue_dictionary.cpp
           )
