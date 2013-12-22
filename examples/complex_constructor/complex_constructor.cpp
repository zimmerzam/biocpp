#include <iostream>
#include <string>
#include <BioCpp.hxx>

#include <BioCpp/pdb/standard/print_atom.hpp>

int main(int argc, char* argv[]){
  if(argc < 2){
    std::cout << "This program shows how to build a simple 'complex' object by using "      
              << "the informations from a pdb file. The program will print all the atoms "  
              << "contained in the complex."                                               
              << std::endl
              << "Usage: ./complex_constructor file.pdb"
              << std::endl
              << "This example is based on BioCpp " 
              << BIOCPP_VERSION_MAJOR << "." << BIOCPP_VERSION_MINOR
              << std::endl
              << std::endl;
    return 1;
  }
  
  const char* filename = argv[1];
  // Initialize a pdb object: store SEQRES, number of models, number of chains per model...
  BioCpp::pdb::file PDB(filename, BioCpp::pdb::INIT_COMPLETE);
  // Read the first model and create an unstructured container of atoms
  BioCpp::standard::base::model all_info = BioCpp::pdb::readModel<BioCpp::standard::base::atom>(PDB, 1);
  // Order atoms in 'all_info' by using informations stored in the pdb: missing residues are 
  //detected and stored as empty residues 
  BioCpp::standard::base::complex_constructor constr(BioCpp::residue::dictionary);
  BioCpp::standard::base::complex cmp = constr(all_info, PDB.RseqRes, PDB.TseqRes);
  
  // Bonus: print all the atoms in the complex by iterating a function that print a single 
  // atom over the whole complex
  BioCpp::Iterate< BioCpp::standard::base::atom >(cmp, BioCpp::pdb::print_atom);
  return 0;
}
