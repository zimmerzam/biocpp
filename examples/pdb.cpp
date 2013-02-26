/*!
    \file pdb.cpp
    \brief Read a pdb file and print some info
    
    In this example we see how to read a pdb file and retrieve some info such as
    total number of models, number of chains, number of residue per chain and
    primary sequence.
*/

#include <iostream>
#include "../src/BioCpp.h"

int main( int argc, char* argv[] ){
  const char* filename = argc>1 ? argv[1] : "2RNM.pdb"; // if a pdb is passed, read that. else read an example pdb
  
  BioCpp::pdb PDB(filename, 0); // read the pdb file. 
  std::cout << filename << std::endl;
  std::cout << "# models: " << PDB.n_models << std::endl;  // number of models
  std::cout << "# chains: " << PDB.RseqRes.size() << std::endl; // chains per model
  std::cout << "**********" << std::endl;
  std::cout << "id  size" << std::endl;
  for( BioCpp::pdb_seqres_record::iterator ch = PDB.RseqRes.begin(); ch != PDB.RseqRes.end(); ++ch ){
    std::cout << ch->first << "    " << ch->second.size() << std::endl; // chain_id and chain size
  }
  std::cout << "**********" << std::endl;
  std::cout << "SEQRES  :" << PDB.TseqRes[PDB.TseqRes.begin()->first] << std::endl; // SEQRES entry for the first chain
  std::cout << "ATOM    :" << PDB.RseqRes[PDB.RseqRes.begin()->first] << std::endl; // readed sequence for the first chain
  std::cout << std::endl;
  return 0;
}
