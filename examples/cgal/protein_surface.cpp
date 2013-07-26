#include <BioCpp.h>
#include <iostream>
#include <list>

struct chain_surface{
  double alpha;
  std::list<Eigen::Vector3d> points, surface;
  
  chain_surface( double radius ){
    alpha=radius;
  }
  
  void operator()( BioCpp::pdb_atom_info atm ){
    points.push_back(atm.coordinate);
  }
  
  void operator()( BioCpp::standard::chain& ch){
    points.clear();  
    BioCpp::Iterate<BioCpp::pdb_atom_info>(ch, *this);
    BioCpp::AlphaExposed( points, surface, alpha );
  }
};

int main( int argc, char* argv[] ){
  
  BioCpp::pdb PDB(argv[1], 0);  
  BioCpp::pdb_model all_info = PDB.getModel(1);
  BioCpp::standard::complex cmp( all_info, PDB.RseqRes, PDB.RseqRes );
  chain_surface surf(4.5);
  BioCpp::Iterate<BioCpp::standard::chain>(cmp,surf);
  
  for(std::list<Eigen::Vector3d>::iterator coord = surf.surface.begin(); coord!=surf.surface.end(); ++coord){
    std::cout << coord->transpose() << std::endl;
  }  
  
  return 0;
}
