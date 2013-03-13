/*!
    \file fasta_align.cpp
    \brief Read two sequences and try to align them
    
    In this example we see how two sequences are aligned.
*/

#include <iostream>
#include <string>
#include "../src/BioCpp.h"

int main(int argc, char* argv[]){
  std::string fasta1 = argc > 1 ? argv[1] : "ABD-G-";
  std::string fasta2 = argc > 2 ? argv[2] : "BDEFGRE";
  fasta2.insert(fasta2.begin(),'-');
  fasta2.insert(fasta2.end(),'-');
  
  BioCpp::error err = BioCpp::ERR_NONE;
  BioCpp::warning war = BioCpp::WAR_NONE;

  double score = BioCpp::fasta::StrictNeedlemanWunsch(fasta1,fasta2, BioCpp::fasta::ZIMM1, err, war);

  std::cout << fasta1 << std::endl 
            << fasta2 << std::endl; // output ABDEFGRE\n  ^BDEFGRE
  std::cout << score << std::endl;
            
  std::cout << "Warnings: " << BioCpp::warning_to_string[war] << std::endl;
  std::cout << "Errors:   " << BioCpp::error_to_string[err] << std::endl;
  return 0;
}
