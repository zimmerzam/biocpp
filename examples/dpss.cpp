/*!
    \file dpss.cpp
    \brief Example usage of dpss utility
    
    The example shows how to determine a protein secondary structure according to dpss definition
*/

#include <BioCpp.h>

// this functor takes two iterators to residue as parameters, check if they 
// belong to a secondary structure and print some info.
struct print_sec{
  BioCpp::standard::h_bridge_map map; // save your h_bridge_map here!

  template <typename T, 
            typename = typename std::enable_if< std::is_same<typename std::iterator_traits<T>::value_type, BioCpp::standard::residue>::value, BioCpp::standard::residue>::type >
  void operator()( T& res1, T& res2 ){
    if( BioCpp::dpss::getSecondaryStructure(map, res1, res2) != BioCpp::dpss::NOT_A_SECONDARY_STRUCTURE ){
				std::cout << (*res1)[BioCpp::atom::CA].chainId << (*res1)[BioCpp::atom::CA].resSeq << "  " 
									<< (*res2)[BioCpp::atom::CA].chainId << (*res2)[BioCpp::atom::CA].resSeq << "  " 
									<< BioCpp::dpss::getSecondaryStructure(map, res1, res2)
									<< std::endl;
	  }
  }
};

int main(){
  BioCpp::pdb::pdb PDB("2RNM.pdb");		// read a pdb file

  BioCpp::pdb::model<BioCpp::pdb::atom_info>::type all_info = PDB.getModel<BioCpp::pdb::atom_info>(1); // get the first model
  BioCpp::standard::complex cmp(all_info, PDB.TseqRes, PDB.TseqRes); // build a complex
  
  print_sec print; // initialize the print_sec functor
  print.map = BioCpp::standard::h_bridge_map(cmp); // compute the h_bridge_map

  BioCpp::Iterate_iter<BioCpp::standard::residue, BioCpp::standard::residue>( cmp, cmp, print ); // apply the functor

  return 0;
}
