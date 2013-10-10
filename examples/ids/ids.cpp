/*!
    \file ids.cpp
    \brief Example usage of element::id, atom::id and amino_acid::id
    
    In this example we see how to use BioCpp identifiers and some variables
    related to them.
*/

#include <iostream>
#include <BioCpp.h>

int main(){

  std::cout << "***********" << std::endl;
  std::cout << "element id: " << std::endl;
  std::cout << BioCpp::element::C << std::endl; // output: " C"
  std::cout << BioCpp::element::string_to_id["HE22"] << std::endl; // output: " H"
  std::cout << BioCpp::element::id_to_string[BioCpp::element::N] << std::endl; // output: " N"
  std::cout << "***********" << std::endl;
  std::cout << "atom id: " << std::endl;
  std::cout << BioCpp::atom::CA << std::endl; // output: " CA "
  std::cout << BioCpp::atom::string_to_id["2HE2"] << std::endl; // output: "HE22"
  std::cout << BioCpp::atom::id_to_string[BioCpp::atom::N_] << std::endl; // output: " N  "
  std::cout << "***********" << std::endl;
  std::cout << "amino acid id: " << std::endl;
  std::cout << BioCpp::amino_acid::ALA << std::endl; // output: "ALA"
  std::cout << BioCpp::amino_acid::string_to_id["GLY"] << std::endl; // output: "GLY"
  std::cout << BioCpp::amino_acid::id_to_string[BioCpp::amino_acid::PRO] << std::endl; // output: "PRO"
  std::cout << BioCpp::amino_acid::id_to_1_letter[BioCpp::amino_acid::PRO] << std::endl; // output: "P"
  std::cout << "***********" << std::endl;
  return 0;
}
