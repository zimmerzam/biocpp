/*!
    \file complex_iterator.cpp
    \brief Read a pdb file and print coordinates and total number of heavy atoms
    
    This program print all the coordiantes of heavy atom in the protein.
    The first line represent the total number of atoms found.
    The output can be read by some cgal example program in order to determine 
    the protein surface
*/

#include <iostream>
#include <BioCpp.hxx>

// This is a very simple functor that can be used to print 
// information about atoms
struct print_coordinate{
  // This function takes a reference to an atom as arguments and, in case
  // it is not an hydrogen, it prints its coordinate
  void operator()(BioCpp::standard::base::atom& atom){
    if(atom.element==BioCpp::element::H){
      return;
    }
    std::cout << atom.coordinate[0] << "  "
              << atom.coordinate[1] << "  "
              << atom.coordinate[2] << std::endl;
  }
};

// This is a more complicated examples of functor that it is used to
// count heavy atoms inside a container
struct count_atoms{
  int atoms;
  // initialize the number of observed atoms to zero
  count_atoms(){
    atoms = 0;
  }
  // each time this operator is called, the number of observed 
  // atoms is increased
  void operator()(BioCpp::standard::base::atom& atom){
    if(atom.element==BioCpp::element::id::H){
      return;
    }
    ++atoms;
  }
};

int main( int argc, char* argv[] ){
  if(argc<2){
    std::cout << "This program shows how to iterate some functions over all the atoms "      
              << "contained in the complex. A similar strategy can be used for iterating "
              << "over all residues (or chains) of a complex. Two kind of functions are"
              << "considered: the first one simply print the coordinates of each atoms;"
              << "the second one update a local variable depending on the kind of atom "
              << "you passed as an argument."
              << std::endl
              << "Usage: ./complex_iterator file.pdb"
              << std::endl
              << "This example is based on BioCpp " 
              << BIOCPP_VERSION_MAJOR << "." << BIOCPP_VERSION_MINOR
              << std::endl
              << std::endl;
    return 1;
  }

  const char* filename = argv[1];
  
  // Initialize a pdb object: store SEQRES, number of models, number of chains per model...
  BioCpp::pdb::file PDB(filename, 0); 
  // Read the first model and create an unstructured container of atoms
  BioCpp::standard::base::model all_info = BioCpp::pdb::readModel<BioCpp::standard::base::atom>(PDB,1);
  // Order atoms in 'all_info' by using informations stored in the pdb: missing residues are 
  //detected and stored as empty residues 
  BioCpp::standard::base::complex_constructor cmp_constr(BioCpp::residue::dictionary);
  BioCpp::standard::base::complex cmp = cmp_constr( all_info, PDB.RseqRes, PDB.RseqRes );
  
  // count all the heavy atoms in the complex
  count_atoms count;
  BioCpp::Iterate<BioCpp::standard::base::atom>(cmp,count);
  std::cout << count.atoms << std::endl;
  
  // iterate a print_coordinate functor over the complex
  print_coordinate printer;
  BioCpp::Iterate<BioCpp::standard::base::atom>(cmp,printer);
  
 return 0;
}
