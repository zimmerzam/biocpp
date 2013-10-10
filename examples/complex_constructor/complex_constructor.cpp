#include <iostream>
#include <string>
#include <BioCpp.h>

int main(int argc, char* argv[]){
  const char* filename = argv[1];
  BioCpp::pdb::pdb PDB(filename, BioCpp::pdb::INIT_COMPLETE);
  BioCpp::standard::base::model all_info = PDB.getModel<BioCpp::standard::base::atom>(1);

  BioCpp::standard::base::complex_constructor constr;

  BioCpp::standard::base::complex cmp = constr(all_info, PDB.RseqRes, PDB.TseqRes);
  
  BioCpp::pdb::print_atom_line printer(std::cout);
  BioCpp::Iterate< BioCpp::standard::base::atom >(cmp, printer);

  return 0;
}
