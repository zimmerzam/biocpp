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
  std::string fasta2 = argc > 2 ? argv[2] : "-BDEFGRE";

  BioCpp::fasta::StrictNeedlemanWunschAlignment align;

  BioCpp::error err = align(fasta1, fasta2);
  
  if(err==0)
    std::cout << fasta1 << std::endl 
              << fasta2 << std::endl; // output ABDEFGRE\n  ^BDEFGRE
  else
    std::cout << "failed" << std::endl;
  return 0;
}
