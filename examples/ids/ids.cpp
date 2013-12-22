/*!
    \file ids.cpp
    \brief Example usage of element::id, atom::id and amino_acid::id
    
    In this example we see how to use BioCpp identifiers and some variables
    related to them.
*/

#include <iostream>
#include <BioCpp.hxx>

int main(){

  std::cout << "***********" << std::endl;
  std::cout << "element id: " << std::endl;
  std::cout << BioCpp::element::C << std::endl; // output: " C"
  std::cout << BioCpp::element::dictionary.string_to_id["HE22"] << std::endl; // output: " H"
  std::cout << BioCpp::element::dictionary.id_to_string[BioCpp::element::N] << std::endl; // output: " N"
  std::cout << "***********" << std::endl;
  std::cout << "atom id: " << std::endl;
  std::cout << BioCpp::atom::CA << std::endl; // output: " CA "
  std::cout << BioCpp::atom::dictionary.string_to_id["2HE2"] << std::endl; // output: "HE22"
  std::cout << BioCpp::atom::dictionary.id_to_string[BioCpp::atom::N] << std::endl; // output: " N  "
  std::cout << "***********" << std::endl;
  std::cout << "amino acid id: " << std::endl;
  std::cout << BioCpp::residue::ALA << std::endl; // output: "ALA"
  std::cout << BioCpp::residue::dictionary.string_to_id["GLY"] << std::endl; // output: "GLY"
  std::cout << BioCpp::residue::dictionary.id_to_string[BioCpp::residue::PRO] << std::endl; // output: "PRO"
  std::cout << BioCpp::residue::dictionary.definition[BioCpp::residue::PRO].one_letter_name << std::endl; // output: "P"
  std::cout << "***********" << std::endl;
  return 0;
}
