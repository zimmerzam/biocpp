set(BIOCPP_LIBRARY_FOLDER ${BIOCPP_LIBRARY_FOLDER} PARENT_SCOPE)
set(LIBRARY_OUTPUT_PATH ${BIOCPP_LIBRARY_FOLDER})

add_library(biocpp-std-proteins base/dpss/h_bridge_map_constructor.cpp
           )
