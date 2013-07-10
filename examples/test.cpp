#include <BioCpp.h>

int main( int argc, char* argv[] ){
  for(int i = 1; i < argc; ++i){
    BioCpp::pdb PDB(argv[i], 0);
  }

  return 0;
}
