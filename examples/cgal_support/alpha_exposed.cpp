#include <BioCpp.h>

int main( int argc, char* argv[] ){
  BioCpp::pdb::pdb PDB(argv[1], 0);  
  double radius = atof(argv[2]);
  BioCpp::standard::morphology::model all_info = PDB.getModel< BioCpp::standard::morphology::atom >(1);
  
  std::list<BioCpp::standard::morphology::alpha_shape_3::Classification_type> classification_type = {BioCpp::standard::morphology::alpha_shape_3::REGULAR, BioCpp::standard::morphology::alpha_shape_3::EXTERIOR, BioCpp::standard::morphology::alpha_shape_3::SINGULAR};
  std::vector<BioCpp::standard::morphology::alpha_shape_3::Vertex_handle> vertices2;
  BioCpp::morphology::alpha_vertices<BioCpp::standard::morphology::alpha_shape_3, BioCpp::standard::morphology::atom>(all_info,radius,vertices2, BioCpp::standard::morphology::alpha_shape_3::GENERAL, classification_type);

  std::size_t nbf2 = vertices2.size();
  for(std::size_t i = 0; i < nbf2; ++i){
    std::cout << vertices2[i]->point().serial << " " 
              << vertices2[i]->point().id << " "
              << vertices2[i]->point().resName << "   "
              << vertices2[i]->point().resSeq << "   "
              << 1
              << std::endl;
  }
  std::cout << "---------------------------------------" << std::endl
            << "Number of surface atoms    =   " << nbf2 << std::endl
            << "Number of buried atoms     =   " << all_info.size()-nbf2 << std::endl;
  return 0;
}
