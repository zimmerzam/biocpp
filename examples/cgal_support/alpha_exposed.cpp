#include <BioCpp.h>

int main( int argc, char* argv[] ){
  BioCpp::pdb::pdb PDB(argv[1], 0);  
  double radius = atof(argv[2]);
  BioCpp::standard::morphology::model all_info = PDB.getModel< BioCpp::standard::morphology::atom >(1);
  
  std::cout << "Find protein surface using cgal functions" << std::endl;
  BioCpp::standard::morphology::alpha_shape_3 as1( all_info.begin(), all_info.end(), radius, BioCpp::standard::morphology::alpha_shape_3::GENERAL );

  std::vector<BioCpp::standard::morphology::alpha_shape_3::Vertex_handle> vertices1;
  as1.get_alpha_shape_vertices( std::back_inserter(vertices1), BioCpp::standard::morphology::alpha_shape_3::REGULAR );
  as1.get_alpha_shape_vertices( std::back_inserter(vertices1), BioCpp::standard::morphology::alpha_shape_3::EXTERIOR );
  as1.get_alpha_shape_vertices( std::back_inserter(vertices1), BioCpp::standard::morphology::alpha_shape_3::SINGULAR );

  std::size_t nbf1 = vertices1.size();
  for(std::size_t i = 0; i < nbf1; ++i){
    std::cout << vertices1[i]->point() << std::endl;
  }
  std::cout << "##########################################" << std::endl;
  std::cout << "Now by using built-in function" << std::endl;
  
  std::list<BioCpp::standard::morphology::alpha_shape_3::Classification_type> classification_type = {BioCpp::standard::morphology::alpha_shape_3::REGULAR, BioCpp::standard::morphology::alpha_shape_3::EXTERIOR, BioCpp::standard::morphology::alpha_shape_3::SINGULAR};
  std::vector<BioCpp::standard::morphology::alpha_shape_3::Vertex_handle> vertices2;
  BioCpp::morphology::alpha_vertices<BioCpp::standard::morphology::alpha_shape_3, BioCpp::standard::morphology::atom>(all_info,radius,vertices2, BioCpp::standard::morphology::alpha_shape_3::GENERAL, classification_type);

  std::size_t nbf2 = vertices2.size();
  for(std::size_t i = 0; i < nbf2; ++i){
    std::cout << vertices2[i]->point() << std::endl;
  }
  return 0;
}
