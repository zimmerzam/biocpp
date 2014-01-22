set(BIOCPP_LIBRARY_FOLDER ${BIOCPP_LIBRARY_FOLDER} PARENT_SCOPE)
set(LIBRARY_OUTPUT_PATH ${BIOCPP_LIBRARY_FOLDER})

add_library(biocpp-proteins dpss/dpss_id.cpp
                            dpss/h_bridge_energy.cpp
                            atom_position_reconstruction/bb_hydrogen.cpp
                            atom_position_reconstruction/bb_oxygen.cpp
                            atom_position_reconstruction/beta_carbon.cpp
           )
