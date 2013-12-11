#include <iostream>
#include <BioCpp.h>

int main( int argc, char* argv[] ){
  if(argc<4){
    std::cout << "This program merge each model in pdb1 with the first one in pdb2 and "      
              << "print the merged pdb file."
              << std::endl
              << "Usage: ./pdb_merge pdb1.pdb pdb2.pdb skip(int)"
              << std::endl
              << "This example is based on BioCpp " 
              << BIOCPP_VERSION_MAJOR << "." << BIOCPP_VERSION_MINOR
              << std::endl
              << std::endl;
    return 1;
  }

  const char* filename1 = argv[1];
  const char* filename2 = argv[2];
  int skip = atoi(argv[3]);
  BioCpp::pdb::pdb PDB1(filename1, 0);
  BioCpp::pdb::pdb PDB2(filename2, 0);
  BioCpp::standard::base::complex_constructor cmp_constr;
  // read the complex to be added to each model in pdb1
  BioCpp::standard::base::model all_info2 = PDB2.getModel<BioCpp::standard::base::atom>(1);
  BioCpp::standard::base::complex cmp2 = cmp_constr( all_info2, PDB2.RseqRes, PDB2.RseqRes );
  // printer
  BioCpp::pdb::print_atom_line printer(std::cout);
  for(int mdl=0; mdl < PDB1.n_models; mdl+=skip){
    BioCpp::standard::base::model all_info = PDB1.getModel<BioCpp::standard::base::atom>(mdl+2);
    BioCpp::standard::base::complex cmp = cmp_constr( all_info, PDB1.RseqRes, PDB1.RseqRes );
    std::cout << "MODEL       " << mdl/skip+1 << std::endl;
    BioCpp::Iterate<BioCpp::standard::base::atom>(cmp,printer);
    BioCpp::Iterate<BioCpp::standard::base::atom>(cmp2,printer);
    std::cout << "ENDMDL" << std::endl;
  }
  
  return 0;
}
