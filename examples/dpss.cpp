/*!
    \file dpss.cpp
    \brief Example usage of dpss utility
    
    The example shows how to determine a protein secondary structure according to dpss definition
*/

#include <BioCpp.h>

// this functor takes two iterators to residue as parameters, check if they 
// belong to a secondary structure and print some info.
struct print_sec{
  BioCpp::standard::base::dpss::h_bridge_map map; // save your h_bridge_map here!

  template <typename T, 
            typename = typename std::enable_if< std::is_same<typename std::iterator_traits<T>::value_type, BioCpp::standard::base::residue>::value, BioCpp::standard::base::residue>::type >
  void operator()( T& res1, T& res2 ){
    if( BioCpp::dpss::getSecondaryStructure(map, res1, res2) != BioCpp::dpss::NOT_A_SECONDARY_STRUCTURE ){
				std::cout << (*res1)[BioCpp::atom::CA].chainId << (*res1)[BioCpp::atom::CA].resSeq << "  " 
									<< (*res2)[BioCpp::atom::CA].chainId << (*res2)[BioCpp::atom::CA].resSeq << "  " 
									<< BioCpp::dpss::getSecondaryStructure(map, res1, res2)
									<< std::endl;
	  }
  }
};

int main(int argc, char* argv[]){
  BioCpp::pdb::pdb PDB(argv[1], BioCpp::pdb::INIT_COMPLETE);		// read a pdb file

  BioCpp::standard::base::model all_info = PDB.getModel<BioCpp::standard::base::atom>(1); // get the first model
  
  BioCpp::standard::base::complex_constructor cmp_constr;
  BioCpp::standard::base::complex cmp = cmp_constr(all_info, PDB.RseqRes, PDB.RseqRes); // build a complex
  
  BioCpp::standard::base::dpss::h_bridge_map_constructor h_map_constr;

  print_sec print; // initialize the print_sec functor
  print.map = h_map_constr(cmp); // compute the h_bridge_map

  BioCpp::Iterate_iter<BioCpp::standard::base::residue, BioCpp::standard::base::residue>( cmp, cmp, print ); // apply the functor

  return 0;
}
