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
  const char* filename = argc>1 ? argv[1] : "2RNM.pdb"; // if a pdb is passed, read that. else read an example pdb
  
  // reading element dictionary...
  libconfig::Config ele_cfg;
  ele_cfg.readFile(argv[2]);
  libconfig::Setting& ele_root = ele_cfg.getRoot();
  element_dictionary eleDict;
  eleDict.importSetting( ele_root);
  // reading atom dictionary...
  libconfig::Config atm_cfg;
  atm_cfg.readFile(argv[3]);
  libconfig::Setting& atm_root = atm_cfg.getRoot();
  atom_dictionary atmDict;
  atmDict.importSetting(atm_root);
  // reading residue dictionary...
  libconfig::Config res_cfg;
  res_cfg.readFile(argv[4]);
  libconfig::Setting& res_root = res_cfg.getRoot();
  residue_dictionary resDict;
  resDict.importSetting(res_root);
  
  BioCpp::io::pdb::file PDB(filename, 0); // read the pdb file. 
	
	model mdl = PDB.readModel<atom>(1);
	print_atom_t printer(std::cout, eleDict, atmDict, resDict);
//	BioCpp::Iterate<atom>(mdl,printer);

  BioCpp::residue::dictionary.writeSetting("out.cfg", "residues");
	
  return 0;
}
