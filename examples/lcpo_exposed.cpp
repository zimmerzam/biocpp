#include <BioCpp.h>

struct exposition{
  double threshold;
  BioCpp::base_neighborood_map<BioCpp::standard::base::atom*> map;
  
  int exposed, buried;
  
  exposition(){
    exposed = 0;
    buried = 0;
  }
  
  void operator()(BioCpp::standard::base::atom& at1, BioCpp::standard::base::atom& at2){
    if(&at1==&at2){
      return;
    }
    double distance = (at1.coordinate-at2.coordinate).norm();
    if( distance < 2.*threshold ){
      map[&at1].insert(&at2);
      map[&at2].insert(&at1);
	  }
  }
  
  void operator()(BioCpp::standard::base::atom& at){
    double area = BioCpp::SurfaceAreaLCPO( at, map, threshold );
    if(area>0){
      ++exposed;
    }
    else{
      ++buried;
    }
    std::cout << at.serial << " " 
              << at.id << " "
              << at.resName << "   "
              << at.resSeq << "   "
              << area
              << std::endl;
  }
  
  void report(){
    std::cout << "---------------------------------------" << std::endl
              << "Number of surface atoms    =   " << exposed << std::endl
              << "Number of buried atoms     =   " << buried << std::endl;
  }
};

int main(int argc, char* argv[]){
  BioCpp::pdb::pdb PDB(argv[1], BioCpp::pdb::INIT_COMPLETE);		// read a pdb file
  BioCpp::standard::base::model all_info = PDB.getModel<BioCpp::standard::base::atom>(1); // get the first model
  BioCpp::standard::base::complex_constructor cmp_constr;
  BioCpp::standard::base::complex cmp = cmp_constr(all_info, PDB.RseqRes, PDB.RseqRes); // build a complex

  double radius = atof(argv[2]);
  
  exposition compute;
  compute.threshold = radius;
  BioCpp::Iterate<BioCpp::standard::base::atom, BioCpp::standard::base::atom>(cmp, cmp, compute);
  BioCpp::Iterate<BioCpp::standard::base::atom>(cmp, compute);
  compute.report();
  
  return 0;
}
