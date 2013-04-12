#include <map>
#include "../src/BioCpp_default.h"

typedef std::pair< char,int > residue;
typedef std::pair< residue, residue > contact;
typedef std::map< contact, unsigned int > contact_map;

struct populate_contacts{
  contact_map contacts;
  bool use_hydrogens;
  double contact_distance;
  void operator()( BioCpp::standard::residue& res1, BioCpp::standard::residue& res2 ){
    for( BioCpp::standard::residue::iterator r1 = res1.begin(); r1 !=res1.end(); ++r1 ){
      for( BioCpp::standard::residue::iterator r2 = res2.begin(); r2 !=res2.end(); ++r2 ){
        if(&res1==&res2){
          continue;
        }
        if( not use_hydrogens){
          if( r1->element==BioCpp::element::H or r2->element==BioCpp::element::H ){
            continue;
          }
        }
        double distance = (r1->coordinate-r2->coordinate).norm();
        if( distance < contact_distance ){
          residue resi1 = std::make_pair( r1->chainId, r1->resSeq );
          residue resi2 = std::make_pair( r2->chainId, r2->resSeq );
          contact cnt = std::make_pair( resi1, resi2 );
          ++contacts[cnt];
          return;
        }
      }
    }
  }
};

std::ostream& operator << (std::ostream& out, contact_map map){
  for( contact_map::iterator it = map.begin(); it != map.end(); ++it ){
    out << it->first.first.first << "  " << it->first.first.second << "  ";
    out << it->first.second.first << "  " << it->first.second.second << "  ";
    out << it->second;
    out << std::endl;
  }
  return out;
}

int main(int argc, char* argv[]){
  const char* filename = argv[1];
  
  populate_contacts populate;
  populate.use_hydrogens = bool(atoi(argv[2]));
  populate.contact_distance = atof(argv[3]);

  BioCpp::pdb PDB(filename, 0);
  for(int mdl = 1; mdl <= PDB.n_models; ++mdl){
    BioCpp::pdb_model all_info = PDB.getModel(mdl);
    BioCpp::standard::complex cmp( all_info, PDB.RseqRes );
    BioCpp::Iterate<BioCpp::standard::residue, BioCpp::standard::residue>(cmp,cmp,populate);
  }

  std::cout << populate.contacts;
  return 0;
}
