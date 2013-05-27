/*!
    \file fasta_align.cpp
    \brief Read two sequences and try to align them
    
    In this example we see how two sequences are aligned.
*/

#include <iostream>
#include <string>
#include <BioCpp.h>

int main(int argc, char* argv[]){
  std::string fasta1 = argc > 1 ? argv[1] : "ABD-G-";
  std::string fasta2 = argc > 2 ? argv[2] : "-BDEFGRE";
  
  BioCpp::error err = BioCpp::ERR_NONE;
  BioCpp::warning war = BioCpp::WAR_NONE;

  BioCpp::fasta::NeedlemanWunsch(fasta1,fasta2, BioCpp::fasta::BLOSUM62, err, war);

  std::cout << fasta1 << std::endl 
            << fasta2 << std::endl; // output ABDEFGRE\n  ^BDEFGRE
            
  std::cout << "Warnings: " << BioCpp::warning_to_string[war] << std::endl;
  std::cout << "Errors:   " << BioCpp::error_to_string[err] << std::endl;
  return 0;
}
