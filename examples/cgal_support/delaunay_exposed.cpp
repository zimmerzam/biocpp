#include <BioCpp.h>

int main( int argc, char* argv[] ){
  BioCpp::pdb::pdb PDB(argv[1], 0);  
  BioCpp::standard::morphology::model all_info = PDB.getModel< BioCpp::standard::morphology::atom >(1);

  std::vector<BioCpp::standard::morphology::triangulation_3::Vertex_handle> vertices2;

  BioCpp::morphology::delaunay_surface<BioCpp::standard::morphology::triangulation_3, BioCpp::standard::morphology::atom>(all_info, vertices2);

  std::size_t nbf2 = vertices2.size();
  for(std::size_t i = 0; i < nbf2; ++i){
    std::cout << vertices2[i]->point().coordinate.transpose() << std::endl;
  }
  return 0;
}


