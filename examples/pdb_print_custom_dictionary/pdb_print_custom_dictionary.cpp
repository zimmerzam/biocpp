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
#include <BioCpp/chemical_component_dictionary/elements/element_dictionary.hpp>

typedef BioCpp::base::atom atom;
typedef BioCpp::io::model<atom>::type model;
typedef BioCpp::element::dictionary_t element_dictionary;
typedef BioCpp::atom::dictionary_t    atom_dictionary;
typedef BioCpp::residue::dictionary_t residue_dictionary;
typedef BioCpp::io::pdb::print_atom_t< element_dictionary, atom_dictionary, residue_dictionary > print_atom_t;

int main( int argc, char* argv[] ){
  const char* pdb_filename = argv[1];
  const char* dic_filename = argv[2];
  const char* moi_filename = argv[3];
  
  // reading dictionaries
  libconfig::Config cfg;
  cfg.readFile(dic_filename);
  libconfig::Setting& root = cfg.getRoot();
  element_dictionary eleDict;
  eleDict.importSetting(root,{"elements"});
  atom_dictionary atmDict;
  atmDict.importSetting(root,{"atoms"});
  residue_dictionary resDict;
  resDict.importSetting(root,{"residues"});
  
  BioCpp::io::pdb::file PDB(pdb_filename, 0); // read the pdb file. 
	
	model mdl = PDB.readModel<atom>(1, eleDict, atmDict, resDict);
	print_atom_t printer(std::cout, eleDict, atmDict, resDict);
//	BioCpp::Iterate<atom>(mdl,printer);

//  BioCpp::residue::dictionary.writeSetting("out.cfg", "residues");
	
  return 0;
}
