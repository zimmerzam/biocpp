/*!
    \file pdb.cpp
    \brief Read a pdb file and print it back
    
    In this example we see how to read a structure from a pdb file as well as how to print
    it back in pdb format
*/

#include <iostream>
#include <BioCpp/base_atom/base_atom.hpp>
#include <BioCpp/standard/fileFormat/pdb/file_reader.hpp>
#include <BioCpp/standard/fileFormat/pdb/printer/print_atom_t.hxx>
#include <BioCpp/base_container/Iterate_single.hxx>

typedef BioCpp::base::atom atom;
typedef BioCpp::io::model<atom>::type model;
typedef BioCpp::element::dictionary_t ele_dict;
typedef BioCpp::atom::dictionary_t atm_dict;
typedef BioCpp::residue::dictionary_t res_dict;
typedef BioCpp::io::pdb::print_atom_t< ele_dict, atm_dict, res_dict > print_atom;

int main( int argc, char* argv[] ){
  const char* filename = argc>1 ? argv[1] : "2RNM.pdb"; // if a pdb is passed, read that. else read an example pdb
  
  BioCpp::io::pdb::file PDB(filename, 0); // read the pdb file. 
  std::cout << filename << std::endl;
	
	model mdl = PDB.readModel<atom>(1);
	print_atom printer(std::cout, BioCpp::element::dictionary, BioCpp::atom::dictionary, BioCpp::residue::dictionary);
	BioCpp::Iterate<atom>(mdl,printer);
	
  return 0;
}
