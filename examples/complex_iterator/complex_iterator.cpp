/*!
    \file print_atom_coordinate.cpp
    \brief Read a pdb file and print coordinates of heavy atoms
    
    This program print all the coordiantes of heavy atom in the protein.
    The first line represent the total number of atoms found.
    The output can be read by alpha_shape.cpp to determine the atoms 
    belonging to the protein surface    
*/

#include <iostream>
#include <BioCpp.h>

struct print_coordinate{
  void operator()(BioCpp::standard::base::atom& atom){
    if(atom.element==BioCpp::element::id::H){
      return;
    }
    std::cout << atom.coordinate[0] << "  "
              << atom.coordinate[1] << "  "
              << atom.coordinate[2] << std::endl;
  }
};

struct count_atoms{
  int atoms;
  count_atoms(){
    atoms = 0;
  }
  void operator()(BioCpp::standard::base::atom& atom){
    if(atom.element==BioCpp::element::id::H){
      return;
    }
    ++atoms;
  }
};

int main( int argc, char* argv[] ){
  const char* filename = argc>1 ? argv[1] : "2RNM.pdb"; // if a pdb is passed, read that. else read an example pdb
  
  BioCpp::pdb::pdb PDB(filename, 0); // read the pdb file. 
  BioCpp::standard::base::model all_info = PDB.getModel<BioCpp::standard::base::atom>(1);
  BioCpp::standard::base::complex_constructor cmp_constr;
  BioCpp::standard::base::complex cmp = cmp_constr( all_info, PDB.RseqRes, PDB.RseqRes );
  count_atoms count;
  BioCpp::Iterate<BioCpp::standard::base::atom>(cmp,count);
  std::cout << count.atoms << std::endl;
  print_coordinate printer;
  BioCpp::Iterate<BioCpp::standard::base::atom>(cmp,printer);
 return 0;
}
