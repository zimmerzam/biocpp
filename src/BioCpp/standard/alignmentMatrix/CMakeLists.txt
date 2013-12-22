set(BIOCPP_LIBRARY_FOLDER ${BIOCPP_LIBRARY_FOLDER} PARENT_SCOPE)
set(LIBRARY_OUTPUT_PATH ${BIOCPP_LIBRARY_FOLDER})

add_library(biocpp-std-align BLOSUM45.cpp
                             BLOSUM62.cpp
                             BLOSUM80.cpp
                             PAM30.cpp
                             PAM70.cpp
                             ZIMM1.cpp
           )
